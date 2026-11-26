"""Benign SV annotations: port of ``AnnotSV-benignsv.tcl``.

* ``checkBenignFiles`` (+ the gnomAD, DGV, DDD, 1000g, dbVar, ClinGen,
  ClinVar, IMH, CMRI, PacBio CoLoRS, FGR and HPRC parsers) -> ``check_benign_files``
* ``checkOverlappedGenesBenignFiles``  -> ``check_overlapped_genes_benign_files``
* ``benignSVannotation``     -> ``BenignSV``    (B_<type>_source/coord/AFmax)
* ``poBenignSVannotation``   -> ``PoBenignSV``  (po_B_<gain|loss>_allG/someG_source/coord)

``BenignAnnotations`` gathers the header and per-line text that
``OrganizeAnnotation`` writes right after the pathogenic SNV/indel columns
(Human only).

Tcl quirks kept on purpose (see tests/python/test_benign_parity.py):

* Every comparison of a field with a number uses Tcl ``expr`` semantics: a
  non-numeric field is compared as a string (so an empty AF is "< 0.01").
* ``checkBenignFiles`` loops over GRCh37, GRCh38 and CHM13 whatever the genome
  build, and then drops the ``B_<type>_*`` columns when the CHM13 sorted BED
  is missing (``$benignDir`` keeps its last loop value).
* 1000g: when bcftools is missing, ``exec`` has already created the empty
  ``.tmp.vcf``, which is then read instead of the download (no SV kept). The
  BED start written is the VCF ``POS`` (not ``POS-1``).
* HPRC and 1000g set one variable per INFO key (``set $key $value``); only the
  keys the procs read back are emulated.
* dbVar: the downloads are deleted only when they gave some lines, and only the
  deletions file triggers the removal of the overlappedGenes files.
* FGR fails (Tcl error, exit code 1) unless the four FGR files exist.
* checkOverlappedGenesBenignFiles fails (Tcl error, exit code 1) right after
  writing an empty overlappedGenes file, when no gene overlaps a benign SV
  type. Its lines come in Tcl hash order (not observable: the file is read
  into an array).
* poBenignSVannotation keeps the genes of the Gain and the Loss
  overlappedGenes files in one array (keyed by the benign SV coordinates).
"""

from __future__ import annotations

import glob
import os
import re
import subprocess
import time
from typing import Dict, List, Optional

from . import tcl
from .config import ConfigExit
from .gene_based import tcl_expr_value
from .general import check_bed, lines_from_file, normalize_svtype, write_text_in_file
from .pathogenic import _ann_dir, _any_col, _delete, _gz_lines, _li, _read_lines, _unix_sort, _wanted

TYPES = ("Loss", "Gain", "Ins", "Inv")


class TclError(ConfigExit):
    """A Tcl runtime error: tclsh prints it and exits with 1."""

    def __init__(self, msg: str):
        super().__init__(msg, 1)


# --------------------------------------------------------------------------
# Tcl expr helpers
# --------------------------------------------------------------------------
def _lt(a, b) -> bool:
    return tcl.compare(a, b) < 0


def _gt(a, b) -> bool:
    return tcl.compare(a, b) > 0


def _num(v):
    n = tcl.to_number(v) if isinstance(v, str) else v
    if n is None:
        raise TclError(f"can't use non-numeric string as operand: \"{v}\"")
    return n


def _f4(v) -> str:
    """``format "%.4f" $v``."""
    n = tcl.to_number(v) if isinstance(v, str) else v
    if n is None:
        raise TclError(f"expected floating-point number but got \"{v}\"")
    return "%.4f" % n


def _minus1(pos: str) -> str:
    """``if {$pos ne 0} {set pos [expr {$pos-1}]}``."""
    return pos if pos == "0" else str(_num(pos) - 1)


def _lsearch_re(items: List[str], pattern: str, start: int = 0) -> int:
    for i in range(max(start, 0), len(items)):
        if re.search(pattern, items[i]):
            return i
    return -1


def _benign_dir(c, build: str) -> str:
    return f"{_ann_dir(c)}/SVincludedInFt/BenignSV/{build}"


def _tmp_file(d: str, svtype: str, build: str) -> str:
    return f"{d}/benign_{svtype}_SV_{build}.tmp.bed"


def _write_all(d: str, build: str, to_write: Dict[str, List[str]], order=TYPES) -> None:
    for t in order:
        if to_write.get(t):
            write_text_in_file("\n".join(to_write[t]), _tmp_file(d, t, build))


def _header_index(path: str, ls: List[str], names) -> Dict[str, int]:
    idx = {}
    for key, name in names:
        if name not in ls:
            raise ConfigExit(f"{path}\nBad header line syntax. {name} column not found - Exit with error", 2)
        idx[key] = ls.index(name)
    return idx


def remove_overlapped_genes_benign_files(c, build: str) -> None:
    """removeOverlappedGenesBenignFiles."""
    d = _benign_dir(c, build)
    for tx in ("RefSeq", "ENSEMBL"):
        for t in TYPES:
            _delete(f"{d}/overlappedGenes_{tx}_in_benign_{t}_SV_{build}.tsv")


# --------------------------------------------------------------------------
# Source parsers (checkXXX_benignFile)
# --------------------------------------------------------------------------
def _check_clinvar(c, build: str) -> None:
    d = _benign_dir(c, build)
    files = glob.glob(f"{d}/clinvar*.vcf.gz")
    if not files:
        return
    downloaded = files[-1]
    remove_overlapped_genes_benign_files(c, build)
    to_write: Dict[str, List[str]] = {t: [] for t in TYPES}
    kind = {"Deletion": "Loss", "Duplication": "Gain", "Insertion": "Ins", "Inversion": "Inv"}
    rx = re.compile(r"ALLELEID=([^;]+)?;.*CLNREVSTAT=([^;]+)?;.*CLNSIG=Benign.*?CLNVC=([^;]+)?")
    min_size = c["SVminSize"]
    for L in _gz_lines(downloaded):
        Ls = L.split("\t")
        if L.startswith("#"):
            continue
        m = rx.search(_li(Ls, 7))
        if not m:
            continue
        alleleid, revstat, clnvc = (x if x is not None else "" for x in m.groups())
        if not re.search("criteria_provided|_multiple_submitters|reviewed_by_expert_panel", revstat):
            continue
        if "no_assertion_criteria_provided" in revstat:
            continue
        length = abs(len(_li(Ls, 3)) - len(_li(Ls, 4)))
        if _lt(length, min_size):
            continue
        chrom = _li(Ls, 0)
        start = _minus1(_li(Ls, 1))
        end = _num(start) + length
        coord = f"{chrom}:{start}-{end}"
        if clnvc in kind:
            to_write[kind[clnvc]].append(f"{chrom}\t{start}\t{end}\tCLN:{alleleid}\t{coord}\tCLN")
    _write_all(d, build, to_write)
    _delete(downloaded)


def _check_clingen(c, build: str) -> None:
    d = _benign_dir(c, build)
    f1 = f"{d}/ClinGen_gene_curation_list_{build}.tsv"
    f2 = f"{d}/ClinGen_region_curation_list_{build}.tsv"
    files = [f for f in (f1, f2) if os.path.exists(f)]
    if not files:
        return
    remove_overlapped_genes_benign_files(c, build)
    to_write: Dict[str, List[str]] = {"Loss": [], "Gain": []}
    for path in files:
        idx: Dict[str, int] = {}
        for L in _read_lines(path):
            Ls = L.split("\t")
            if re.search("(^#Gene Symbol)|(^#ISCA ID)", L):
                idx = {"id": 0}
                for key, name, label in (("coord", "Genomic Location", "Genomic Location"),
                                         ("hi", "Haploinsufficiency Score", "Haploinsufficiency Score"),
                                         ("ts", "Triplosensitivity Score", "Triplosensitivity")):
                    if name not in Ls:
                        raise ConfigExit(f"{path}\nBad header line syntax. {label} column not found - Exit with error", 2)
                    idx[key] = Ls.index(name)
                continue
            if L.startswith("#"):
                continue
            if not idx:
                raise TclError("can't read \"i_id\": no such variable")
            ID = _li(Ls, idx["id"])
            coord = re.sub("(chr)| ", "", _li(Ls, idx["coord"]), count=1)
            m = re.search(r"([0-9XYMT]+):([0-9]+)-([0-9]+)", coord)
            if not m:
                continue
            chrom, start, end = m.groups()
            if _li(Ls, idx["hi"]) == "40":
                to_write["Loss"].append(f"{chrom}\t{start}\t{end}\tHI40:{ID}\t{coord}\tHI")
            if _li(Ls, idx["ts"]) == "40":
                to_write["Gain"].append(f"{chrom}\t{start}\t{end}\tTS40:{ID}\t{coord}\tTS")
    _write_all(d, build, to_write)
    _delete(f1)
    _delete(f2)


def _check_dgv(c, build: str) -> None:
    d = _benign_dir(c, build)
    files = glob.glob(f"{d}/GRCh3*_hg*_variants_*.txt")
    if not files:
        return
    path = files[0]
    remove_overlapped_genes_benign_files(c, build)
    to_write: Dict[str, List[str]] = {t: [] for t in TYPES}
    idx: Dict[str, int] = {}
    for L in lines_from_file(path):
        Ls = L.split("\t")
        if L.startswith("variantaccession"):
            idx = {"id": 0}
            idx.update(_header_index(path, Ls, (("chr", "chr"), ("start", "start"), ("end", "end"),
                                                ("n", "samplesize"), ("gains", "observedgains"),
                                                ("loss", "observedlosses"))))
            continue
        if not idx:
            raise TclError("can't read \"i_samplesize\": no such variable")
        samplesize = _li(Ls, idx["n"])
        if not _gt(samplesize, c["minTotalNumber"]):
            continue

        def freq(obs):
            a, b = _num(obs), _num(samplesize)
            a = a * 100
            b = 2 * b
            if isinstance(a, int) and isinstance(b, int):
                if b == 0:
                    raise TclError("divide by zero")
                return a // b
            return a / b

        freqgain, freqloss = freq(_li(Ls, idx["gains"])), freq(_li(Ls, idx["loss"]))
        chrom, start, end = _li(Ls, idx["chr"]), _li(Ls, idx["start"]), _li(Ls, idx["end"])
        head = f"{chrom}\t{start}\t{end}\t{_li(Ls, idx['id'])}\t{chrom}:{start}-{end}\t"
        if freqgain > 0.1:
            to_write["Gain"].append(head + _f4(freqgain * 1.0 / 100))
        if freqloss > 0.1:
            to_write["Loss"].append(head + _f4(freqloss * 1.0 / 100))
    _write_all(d, build, to_write)
    _delete(path)


def _check_gnomad(c, build: str) -> None:
    d = _benign_dir(c, build)
    to_write: Dict[str, List[str]] = {t: [] for t in TYPES}
    if build == "GRCh37":
        path = f"{d}/gnomad_v2.1_sv.sites.bed.gz"
        if not os.path.exists(path):
            return
        remove_overlapped_genes_benign_files(c, build)
        idx: Dict[str, int] = {}
        for L in _gz_lines(path):
            if L == "":
                continue
            Ls = L.split("\t")
            if L.startswith("#chrom"):
                idx = _header_index(path, Ls, (("chrom", "#chrom"), ("start", "start"), ("end", "end"),
                                               ("svid", "name"), ("svtype", "svtype"), ("nhomalt", "N_HOMALT"),
                                               ("popmaxaf", "POPMAX_AF"), ("filter", "FILTER")))
                continue
            if not idx:
                raise TclError("can't read \"i_filter\": no such variable")
            if _li(Ls, idx["filter"]) != "PASS":
                continue
            svtype = re.sub(":.+", "", _li(Ls, idx["svtype"]), count=1)
            if svtype not in ("DEL", "DUP", "INS", "INV"):
                continue
            nhomalt, popmaxaf = _li(Ls, idx["nhomalt"]), _li(Ls, idx["popmaxaf"])
            if _lt(popmaxaf, 0.001) or _lt(nhomalt, 5):
                continue
            chrom, start, end = _li(Ls, idx["chrom"]), _li(Ls, idx["start"]), _li(Ls, idx["end"])
            infos = f"{chrom}\t{start}\t{end}\t{_li(Ls, idx['svid'])}\t{chrom}:{start}-{end}\t{_f4(popmaxaf)}"
            to_write[{"DEL": "Loss", "DUP": "Gain", "INS": "Ins", "INV": "Inv"}[svtype]].append(infos)
        _write_all(d, build, to_write)
        _delete(path)
    elif build == "GRCh38":
        files = glob.glob(f"{d}/gnomad.v4.0.sv.chr*.vcf.gz")
        if not files:
            return
        remove_overlapped_genes_benign_files(c, build)
        for path in files:
            for L in _gz_lines(path):
                if L[:1] == "#" or L == "":
                    continue
                Ls = L.split("\t")
                if _li(Ls, 6) != "PASS":
                    continue
                chrom = _li(Ls, 0)
                pos = _minus1(_li(Ls, 1))
                svid = _li(Ls, 2)
                infos_l = _li(Ls, 7).split(";")
                vals = []
                k = 0
                for key in ("END", "SVTYPE", "controls_and_biobanks_AN", "controls_and_biobanks_AF",
                            "controls_and_biobanks_N_HOMALT"):
                    k = _lsearch_re(infos_l, f"^{key}=", k)
                    if k == -1:
                        break
                    vals.append(_li(infos_l[k].split("="), 1))
                    if key == "SVTYPE" and vals[-1] in ("BND", "CPX", "CTX"):
                        break
                    if key == "controls_and_biobanks_AN" and _lt(vals[-1], 1000):
                        break
                    if key == "controls_and_biobanks_AF" and _lt(vals[-1], 0.001):
                        break
                    if key == "controls_and_biobanks_N_HOMALT" and _lt(vals[-1], 5):
                        break
                else:
                    end, svtype, _, af, _ = vals
                    infos = f"{chrom}\t{pos}\t{end}\t{svid}\t{chrom}:{pos}-{end}\t{_f4(af)}"
                    if "DEL" in svtype:
                        to_write["Loss"].append(infos)
                    if "DUP" in svtype:
                        to_write["Gain"].append(infos)
                    if "INS" in svtype:
                        to_write["Ins"].append(infos)
                    if svtype == "INV":
                        to_write["Inv"].append(infos)
        _write_all(d, build, to_write)
        for path in files:
            _delete(path)


def _check_hprc(c, build: str) -> None:
    if build != "GRCh38":
        return
    d = _benign_dir(c, build)
    path = f"{d}/hprc.GRCh38.pbsv.vcf.gz"
    if not os.path.exists(path):
        return
    remove_overlapped_genes_benign_files(c, build)
    to_write: Dict[str, List[str]] = {t: [] for t in TYPES}
    for L in _gz_lines(path):
        if L[:1] == "#" or L == "":
            continue
        Ls = L.split("\t")
        if _li(Ls, 6) != "PASS":
            continue
        chrom = _li(Ls, 0)
        pos = _minus1(_li(Ls, 1))
        svid = _li(Ls, 2)
        var: Dict[str, str] = {}
        for inf in _li(Ls, 7).split(";"):
            kv = inf.split("=")
            var[kv[0]] = _li(kv, 1)
        if not all(k in var for k in ("END", "SVTYPE", "AF")):
            continue
        svtype, af = var["SVTYPE"], var["AF"]
        if svtype == "BND":
            continue
        if _lt(af, 0.05):
            continue
        infos = f"{chrom}\t{pos}\t{var['END']}\tHPRC:{svid}\t{chrom}:{pos}-{var['END']}\t{_f4(af)}"
        if "DUP" in svtype:
            to_write["Gain"].append(infos)
        if "DEL" in svtype:
            to_write["Loss"].append(infos)
        if "INS" in svtype:
            to_write["Ins"].append(infos)
        if svtype == "INV":
            to_write["Inv"].append(infos)
    _write_all(d, build, to_write, ("Gain", "Loss", "Ins", "Inv"))
    _delete(path)


def _check_ddd(c, build: str) -> None:
    d = _benign_dir(c, build)
    path = f"{d}/population_cnv_{build.lower()}.txt.gz"
    if not os.path.exists(path):
        return
    remove_overlapped_genes_benign_files(c, build)
    to_write: Dict[str, List[str]] = {"Loss": [], "Gain": []}
    idx: Dict[str, int] = {}
    for L in _gz_lines(path):
        Ls = L.split("\t")
        if L.startswith("#population_cnv_id"):
            idx = {"id": 0}
            for key, name in (("chr", "chr"), ("start", "start"), ("end", "end"), ("del", "deletion_frequency"),
                              ("dup", "duplication_frequency"), ("n", "sample_size")):
                idx[key] = _lsearch_re(Ls, name)
                if idx[key] == -1:
                    raise ConfigExit(f"Bad syntax into {path}.\n{name} field not found - Exit with error", 2)
            continue
        if not idx:
            raise TclError("can't read \"i_sample_size\": no such variable")
        if _lt(_li(Ls, idx["n"]), c["minTotalNumber"]):
            continue
        ID = f"DDD:{_li(Ls, idx['id'])}"
        chrom, start, end = _li(Ls, idx["chr"]), _li(Ls, idx["start"]), _li(Ls, idx["end"])
        delfreq, dupfreq = _li(Ls, idx["del"]), _li(Ls, idx["dup"])
        head = f"{chrom}\t{start}\t{end}\t{ID}\t{chrom}:{start}-{end}\t"
        if _gt(delfreq, 0.001):
            to_write["Loss"].append(head + _f4(delfreq))
        if _gt(dupfreq, 0.001):
            to_write["Gain"].append(head + _f4(dupfreq))
    _write_all(d, build, to_write)
    _delete(path)


def _bcftools_norm(c, src: str, dst: str) -> str:
    """``catch {eval exec $bcftools norm -m -both $src > $dst} Message``: the message."""
    with open(dst, "wb") as out:
        try:
            r = subprocess.run([c.get("bcftools", "bcftools"), "norm", "-m", "-both", src],
                               stdout=out, stderr=subprocess.PIPE)
        except OSError:
            return f"couldn't execute \"{c.get('bcftools', 'bcftools')}\": no such file or directory"
    err = r.stderr.decode(errors="replace")
    err = err[:-1] if err.endswith("\n") else err
    if err:
        return err
    return "child process exited abnormally" if r.returncode else ""


def _check_1000g(c, build: str) -> None:
    d = _benign_dir(c, build)
    files = glob.glob(f"{d}/ALL.wgs.mergedSV*.vcf.gz")
    if not files:
        return
    downloaded = files[0]
    remove_overlapped_genes_benign_files(c, build)
    to_write: Dict[str, List[str]] = {t: [] for t in TYPES}
    tmp = re.sub(r".vcf.gz", ".tmp.vcf", downloaded, count=1)
    msg = _bcftools_norm(c, downloaded, tmp)
    if re.search("error|fail", msg, re.I):
        _delete(tmp)
    lines = _read_lines(tmp) if os.path.exists(tmp) else _gz_lines(downloaded)
    min_size = c["SVminSize"]
    for L in lines:
        if L[:1] == "#" or L == "":
            continue
        Ls = L.split("\t")
        if _li(Ls, 6) != "PASS":
            continue
        chrom, pos = _li(Ls, 0), _li(Ls, 1)
        if pos != "0":
            _num(pos)
        ref, alt = _li(Ls, 3), _li(Ls, 4)
        if not re.search(r"<|\[|\]", alt):
            if _lt(abs(len(ref) - len(alt)), min_size):
                continue
        mx = -1
        var: Dict[str, str] = {"SVTYPE": "CNV"}
        for inf in _li(Ls, 7).split(";"):
            kv = inf.split("=")
            key, val = kv[0], _li(kv, 1)
            var[key] = val
            if re.search("AF$", key) and _gt(val, mx):
                mx = val
        if _lt(mx, 0.001):
            continue
        end = var["END"] if "END" in var else str(_num(pos) + 1)
        svtype = var["SVTYPE"]
        if svtype == "CNV":
            svtype = alt
        if re.search("^<CN", alt):
            svtype = alt
        infos = f"{chrom}\t{pos}\t{end}\t1000g\t{chrom}:{pos}-{end}\t{_f4(mx)}"
        if re.search("<CN0>|DEL", svtype):
            to_write["Loss"].append(infos)
        if re.search("<CN[1-9]|ALU|LINE1|SVA", svtype):
            to_write["Gain"].append(infos)
        if svtype == "INS":
            to_write["Ins"].append(infos)
        if svtype == "INV":
            to_write["Inv"].append(infos)
    _write_all(d, build, to_write)
    _delete(downloaded)
    _delete(tmp)


def _check_dbvar(c, build: str) -> None:
    d = _benign_dir(c, build)
    to_write: Dict[str, List[str]] = {"Loss": [], "Gain": [], "Ins": []}
    paths = {}
    for t, kind in (("Loss", "deletions"), ("Gain", "duplications"), ("Ins", "insertions")):
        path = paths[t] = f"{d}/{build}.nr_{kind}.common.bed.gz"
        if not os.path.exists(path):
            continue
        if t == "Loss":
            remove_overlapped_genes_benign_files(c, build)
        for L in _gz_lines(path):
            if L == "":
                continue
            Ls = L.split("\t")
            chrom, start, end = _li(Ls, 0), _li(Ls, 1), _li(Ls, 2)
            to_write[t].append(f"{chrom}\t{start}\t{end}\tdbVar\t{chrom}:{start}-{end}\t0.01")
    for t in ("Loss", "Gain", "Ins"):
        if to_write[t]:
            write_text_in_file("\n".join(to_write[t]), _tmp_file(d, t, build))
            _delete(paths[t])


def _check_imh(c, build: str) -> None:
    d = _benign_dir(c, build)
    files = glob.glob(f"{d}/*callset.public.bedpe.gz")
    if not files:
        return
    path = files[0]
    remove_overlapped_genes_benign_files(c, build)
    to_write: Dict[str, List[str]] = {t: [] for t in TYPES}
    idx: Dict[str, int] = {}
    for L in _gz_lines(path):
        Ls = L.split("\t")
        if L.startswith("#CHROM_A"):
            idx = {"chrA": 0}
            for key, name in (("chrB", "CHROM_B"), ("start", "START_A"), ("end", "END_B"), ("svtype", "TYPE"),
                              ("ann", "INFO_A")):
                idx[key] = _lsearch_re(Ls, name)
                if idx[key] == -1:
                    raise ConfigExit(f"Bad syntax into {path}.\n{name} field not found - Exit with error", 2)
        if L[:1] == "#" or L == "":
            continue
        if not idx:
            raise TclError("can't read \"i_chrA\": no such variable")
        chrom_a = _li(Ls, idx["chrA"]).replace("chr", "", 1)
        chrom_b = _li(Ls, idx["chrB"]).replace("chr", "", 1)
        if chrom_a != chrom_b:
            continue
        m = re.search(";AF=(.+?);", _li(Ls, idx["ann"]))
        if not m:
            continue
        af = m.group(1)
        if _lt(af, 0.001):
            continue
        start, end = _li(Ls, idx["start"]), _li(Ls, idx["end"])
        svtype = normalize_svtype(_li(Ls, idx["svtype"]))
        infos = f"{chrom_a}\t{start}\t{end}\tIMH\t{chrom_a}:{start}-{end}\t{_f4(af)}"
        t = {"DEL": "Loss", "DUP": "Gain", "INS": "Ins", "INV": "Inv"}.get(svtype)
        if t:
            to_write[t].append(infos)
    _write_all(d, build, to_write)
    _delete(path)
    # Tcl: relative to the working directory
    for pattern in ("Build*.public.v2.vcf.gz", "Supplementary_File_*.zip"):
        found = glob.glob(pattern)
        if len(found) == 1:
            _delete(found[0])


def _vcf_freq_source(c, build: str, path: str, rx, groups, min_freq, label, types_kept) -> None:
    """Shared body of checkCMRI_benignFile and checkPacBioCoLoRS_benignFile."""
    d = _benign_dir(c, build)
    remove_overlapped_genes_benign_files(c, build)
    to_write: Dict[str, List[str]] = {t: [] for t in TYPES}
    lines = _gz_lines(path) if path.endswith(".gz") else _read_lines(path)
    for L in lines:
        Ls = L.split("\t")
        if L.startswith("#"):
            continue
        m = rx.search(_li(Ls, 7))
        if not m:
            continue
        g = dict(zip(groups, (x if x is not None else "" for x in m.groups())))
        if _lt(g["total"], 500):
            continue
        if _lt(abs(_num(g["svlen"])), c["SVminSize"]):
            continue
        if _lt(g["freq"], min_freq):
            continue
        chrom = _li(Ls, 0).replace("chr", "", 1)
        start = _minus1(_li(Ls, 1))
        end = g["end"]
        coord = f"{chrom}:{start}-{end}"
        if end == start:
            end = str(_num(start) + 1)
        if _lt(end, start):
            start, end = end, start
        name = label(Ls)
        t = {"DEL": "Loss", "DUP": "Gain", "INS": "Ins", "INV": "Inv"}.get(g["svtype"])
        if t in types_kept:
            to_write[t].append(f"{chrom}\t{start}\t{end}\t{name}\t{coord}\t{_f4(g['freq'])}")
    _write_all(d, build, to_write)
    _delete(path)


def _check_cmri(c, build: str) -> None:
    files = glob.glob(f"{_benign_dir(c, build)}/pb_joint_merged.sv.vcf")
    if not files:
        return
    rx = re.compile(r"END=([^;]+)?;.*SVTYPE=([^;]+)?;.*SVLEN=([^;]+)?;.*SVN=([^;]+)?;.*SVF=([^;]+)?")
    _vcf_freq_source(c, build, files[-1], rx, ("end", "svtype", "svlen", "total", "freq"), 0.001,
                     lambda Ls: f"CMRI:{_li(Ls, 2)}", TYPES)


def _check_pacbio_colors(c, build: str) -> None:
    files = glob.glob(f"{_benign_dir(c, build)}/CoLoRSdb.{build}.v1.0.0.pbsv.jasmine.vcf.gz")
    if not files:
        return
    rx = re.compile(r"SVTYPE=([^;]+)?;END=([^;]+)?;.*SVLEN=([^;]+)?;.*AN=([^;]+)?;.*AF=([^;]+)?")
    _vcf_freq_source(c, build, files[-1], rx, ("svtype", "end", "svlen", "total", "freq"), 0.05,
                     lambda Ls: "PacBioCoLoRS", ("Loss", "Ins", "Inv"))


def _check_fgr(c, build: str) -> None:
    d = _benign_dir(c, build)
    files = glob.glob(f"{d}/FGR_benign_*_SV_{build}.bed")
    if not files:
        return
    remove_overlapped_genes_benign_files(c, build)
    missing = None
    for t in ("Gain", "Loss", "Ins", "Inv"):
        src = f"{d}/FGR_benign_{t}_SV_{build}.bed"
        if os.path.exists(src):
            write_text_in_file("\n".join(lines_from_file(src)), _tmp_file(d, t, build))
        elif missing is None:
            missing = t
    if missing is not None:
        # Tcl: puts "($nLoss SV Loss + $nGain SV Gain)" with an unset variable
        first = next(t for t in ("Loss", "Gain", "Ins", "Inv")
                     if not os.path.exists(f"{d}/FGR_benign_{t}_SV_{build}.bed"))
        raise TclError(f"can't read \"n{first}\": no such variable")
    for f in files:
        _delete(f)


# --------------------------------------------------------------------------
# checkBenignFiles / checkOverlappedGenesBenignFiles
# --------------------------------------------------------------------------
def check_benign_files(config) -> None:
    """checkBenignFiles: build ``benign_<type>_SV_<build>.sorted.bed``."""
    c = config
    d = ""
    for build in ("GRCh37", "GRCh38", "CHM13"):
        d = _benign_dir(c, build)
        for t in TYPES:
            _delete(_tmp_file(d, t, build))
        for check in (_check_gnomad, _check_dgv, _check_ddd, _check_1000g, _check_dbvar, _check_clingen,
                      _check_clinvar, _check_imh, _check_cmri, _check_pacbio_colors, _check_fgr, _check_hprc):
            check(c, build)
        for t in TYPES:
            sorted_bed = f"{d}/benign_{t}_SV_{build}.sorted.bed"
            tmp = _tmp_file(d, t, build)
            if os.path.exists(tmp):
                if os.path.exists(sorted_bed):
                    with open(sorted_bed, "rb") as src:
                        lines = src.read().split(b"\n")
                    if lines[-1] == b"":
                        lines.pop()
                    with open(tmp, "ab") as dst:
                        dst.write(b"".join(line + b"\n" for line in lines))
                try:  # catch {checkBed ...}
                    check_bed(tmp, c["outputDir"], d)
                except ConfigExit:
                    raise
                except Exception:
                    pass
                _delete(tmp)
        for t in TYPES:
            formatted = f"{d}/benign_{t}_SV_{build}.tmp.formatted.bed"
            if not os.path.exists(formatted):
                continue
            _unix_sort(formatted, f"{d}/benign_{t}_SV_{build}.sorted.bed", "checkBenignFiles")
            _delete(formatted)
        for t in ("Gain", "Loss", "Ins", "Inv"):
            _delete(f"{d}/overlappedGenes_in_benign_{t}_SV_{build}.tsv")

    # Tcl: $benignDir is the one of the last genome build (CHM13)
    for t in ("Gain", "Loss", "Ins", "Inv"):
        if not os.path.exists(f"{d}/benign_{t}_SV_CHM13.sorted.bed"):
            c["outputColHeader"] = [e for e in _wanted(c) if not re.search(f"^B_{t.lower()}_", e)]


def _bedtools_intersect(c, a: str, b: str, out_file: str, opts: List[str], caller: str) -> None:
    """``exec bedtools intersect -sorted -a $a -b $b <opts> > $out``, then without ``-sorted``."""
    bt = c.get("bedtools", "bedtools")

    def run(sorted_opt):
        with open(out_file, "w") as out:
            r = subprocess.run([bt, "intersect"] + sorted_opt + ["-a", a, "-b", b] + opts,
                               stdout=out, stderr=subprocess.PIPE)
        return r.returncode == 0 and not r.stderr, r.stderr.decode(errors="replace").rstrip("\n")

    ok, _ = run(["-sorted"])
    if ok:
        return
    ok, err = run([])
    if not ok:
        raise ConfigExit(f"-- {caller} --\n{bt} intersect -sorted -a {a} -b {b} {' '.join(opts)} > {out_file}\n"
                         f"{err or 'child process exited abnormally'}\nExit with error", 2)


def check_overlapped_genes_benign_files(config) -> None:
    """checkOverlappedGenesBenignFiles: ``overlappedGenes_<tx>_in_benign_<type>_SV_<build>.tsv``."""
    c = config
    if c["organism"] != "Human":
        return
    build = c["genomeBuild"]
    d = _benign_dir(c, build)
    genes_dir = f"{_ann_dir(c)}/Genes/{build}"
    for tx in ("RefSeq", "ENSEMBL"):
        genes_bed = f"{genes_dir}/genes.{tx}.sorted.bed"
        genes_bed = genes_bed if os.path.exists(genes_bed) else ""
        for t in TYPES:
            benign_bed = f"{d}/benign_{t}_SV_{build}.sorted.bed"
            benign_bed = benign_bed if os.path.exists(benign_bed) else ""
            out = f"{d}/overlappedGenes_{tx}_in_benign_{t}_SV_{build}.tsv"
            if os.path.exists(out):
                continue
            print(f"\t   >>> creation of {out}")
            tmp = f"{c['outputDir']}/{int(time.time())}_overlappedGenes_{tx}_in_benign_{t}_SV_{build}.tmp.bed"
            _delete(tmp)
            _bedtools_intersect(c, genes_bed, benign_bed, tmp, ["-wa", "-wb"],
                                f"checkOverlappedGenesBenignFiles, {build}, {tx}, {t}")
            genes: Dict[str, List[str]] = {}
            for L in _read_lines(tmp):
                if L == "":
                    continue
                Ls = L.split("\t")
                genes.setdefault(Ls[-2] if len(Ls) > 1 else "", []).append(_li(Ls, 4))
            _delete(tmp)
            # Tcl writes the coordinates in hash order (the file is only read back into an array)
            write_text_in_file("\n".join(f"{coord}\t{';'.join(sorted(set(g)))}" for coord, g in genes.items()), out)
            if not genes:
                raise TclError("can't unset \"L_genes\": no such variable")


# --------------------------------------------------------------------------
# benignSVannotation / poBenignSVannotation
# --------------------------------------------------------------------------
def _out_tmp(c, repl: str) -> str:
    name = re.sub(r"(.formatted)?.bed$", lambda _m: repl, c["bedFile"], count=1, flags=re.I)
    return f"{c['outputDir']}/{os.path.basename(name)}"


def _skip_af(af: str, benign_af) -> bool:
    """``[string is double $benign_AF] && $benign_AF < $g_AnnotSV(benignAF)``."""
    return tcl.string_is_double(af) and _lt(af, benign_af)


class BenignSV:
    """benignSVannotation: B_<type>_source, _coord, _AFmax (benign SV covering the SV)."""

    def __init__(self, config):
        self.config = config
        self.text: Optional[Dict[str, str]] = None
        self.empty = ""

    def _load(self) -> None:
        c = self.config
        build = c["genomeBuild"]
        d = _benign_dir(c, build)
        types = [t for t in ("Gain", "Loss", "Ins", "Inv") if _any_col(c, f"^B_{t.lower()}_")]
        self.empty = "\t".join("" for _ in types for _ in range(3))
        all_sv: List[str] = []
        coord: Dict[tuple, List[str]] = {}
        source: Dict[tuple, List[str]] = {}
        afmax: Dict[tuple, str] = {}
        for t in types:
            bed = f"{d}/benign_{t}_SV_{build}.sorted.bed"
            bed = bed if os.path.exists(bed) else ""
            tmp = _out_tmp(c, f".intersect.benign-{t}")
            _delete(tmp)
            _bedtools_intersect(c, c["fullAndSplitBedFile"], bed, tmp, ["-f", "1.0", "-wa", "-wb"],
                                f"benignSVAnnotation, {t}")
            for L in _read_lines(tmp):
                if L == "":
                    continue
                Ls = L.split("\t")
                af = Ls[-1]
                if _skip_af(af, c["benignAF"]):
                    continue
                key = f"{_li(Ls, 0)},{_li(Ls, 1)},{_li(Ls, 2)}"
                k = (key, t)
                all_sv.append(key)
                coord.setdefault(k, []).append(_li(Ls, len(Ls) - 2))
                source.setdefault(k, []).append(_li(Ls, len(Ls) - 3))
                if k not in afmax:
                    afmax[k] = af if tcl.string_is_double(af) else ""
                elif tcl.string_is_double(af) and _gt(af, afmax[k]):
                    afmax[k] = af
            _delete(tmp)
        fr = c["metrics"] == "fr"
        self.text = {}
        for key in sorted(set(all_sv)):
            parts: List[str] = []
            for t in types:
                k = (key, t)
                if k in coord:
                    af = afmax[k]
                    if fr:
                        af = tcl_expr_value(af).replace(".", ",")
                    parts += [";".join(sorted(set(source[k]))), ";".join(coord[k]), af]
                else:
                    parts += ["", "", ""]
            self.text[key] = "\t".join(parts)

    def __call__(self, chrom: str, start: str, end: str) -> str:
        if self.text is None:
            self._load()
        return self.text.get(f"{chrom},{start},{end}", self.empty)


class PoBenignSV:
    """poBenignSVannotation: po_B_<gain|loss>_allG_source/_coord, po_B_<gain|loss>_someG_source/_coord."""

    def __init__(self, config):
        self.config = config
        self.coord_source: Optional[Dict[tuple, List[str]]] = None
        self.genes: Dict[str, List[str]] = {}

    def _types(self) -> List[str]:
        return [t for t in ("Gain", "Loss") if _any_col(self.config, f"^po_B_{t.lower()}_")]

    def _load(self) -> None:
        c = self.config
        build = c["genomeBuild"]
        d = _benign_dir(c, build)
        self.coord_source = {}
        for t in self._types():
            ov = f"{d}/overlappedGenes_{c['tx']}_in_benign_{t}_SV_{build}.tsv"
            if not os.path.exists(ov):
                raise TclError(f"couldn't open \"{ov}\": no such file or directory")
            for L in lines_from_file(ov):
                self.genes.setdefault(tcl.lindex(L, 0), []).append(tcl.lindex(L, 1))
            bed = f"{d}/benign_{t}_SV_{build}.sorted.bed"
            bed = bed if os.path.exists(bed) else ""
            tmp = _out_tmp(c, f".intersect.po_benign-{t}")
            _delete(tmp)
            _bedtools_intersect(c, c["bedFile"], bed, tmp, ["-wa", "-wb"], f"benignSVAnnotation, {t}")
            for L in _read_lines(tmp):
                if L == "":
                    continue
                Ls = L.split("\t")
                if _skip_af(Ls[-1], c["benignAF"]):
                    continue
                key = (f"{_li(Ls, 0)},{_li(Ls, 1)},{_li(Ls, 2)}", t)
                self.coord_source.setdefault(key, []).append(f"{_li(Ls, len(Ls) - 2)} {_li(Ls, len(Ls) - 3)}")
            # Tcl deletes the intersection at the end of OrganizeAnnotation
            _delete(tmp)

    def __call__(self, chrom: str, start: str, end: str, genes_sv: str) -> str:
        if self.coord_source is None:
            self._load()
        sv = f"{chrom},{start},{end}"
        parts: List[str] = []
        for t in self._types():
            groups: Dict[str, List[tuple]] = {"allG": [], "someG": []}
            for cs in self.coord_source.get((sv, t), []):
                items = tcl.split_list(cs)
                b_coord, b_source = _li(items, 0), _li(items, 1)
                which = "someG"
                if b_coord in self.genes and genes_sv != "":
                    benign_genes = ";".join(self.genes[b_coord]).split(";")
                    if all(g in benign_genes for g in genes_sv.split(";")):
                        which = "allG"
                groups[which].append((b_source, b_coord))
            for which in ("allG", "someG"):
                g = groups[which]
                if not g:
                    parts += ["", ""]
                elif len(g) > 20:
                    parts += [";".join(s for s, _ in g[:20]) + "...", ";".join(x for _, x in g[:20]) + "..."]
                else:
                    parts += [";".join(s for s, _ in g), ";".join(x for _, x in g)]
        return "\t".join(parts)


# --------------------------------------------------------------------------
# OrganizeAnnotation part
# --------------------------------------------------------------------------
class BenignAnnotations:
    """The benign columns that OrganizeAnnotation writes after the pathogenic SNV/indel ones (Human only)."""

    def __init__(self, config):
        self.config = config
        self.sv = BenignSV(config)
        self.po = PoBenignSV(config)

    @property
    def human(self) -> bool:
        return self.config["organism"] == "Human"

    def header(self) -> str:
        c = self.config
        if not self.human:
            return ""
        h = ""
        for t in ("gain", "loss", "ins", "inv"):
            if _any_col(c, f"^B_{t}_"):
                h += f"\tB_{t}_source\tB_{t}_coord\tB_{t}_AFmax"
        for t in ("gain", "loss"):
            if _any_col(c, f"^po_B_{t}_"):
                h += f"\tpo_B_{t}_allG_source\tpo_B_{t}_allG_coord\tpo_B_{t}_someG_source\tpo_B_{t}_someG_coord"
        return h

    def text(self, chrom: str, start: str, end: str, sv_left: str, sv_right: str, mode: str,
             gene_name: str) -> str:
        """``start``/``end``: the intersection on split lines, the SV on full lines."""
        if not self.human:
            return ""
        if mode == "split":
            po = self.po("NULL", "NULL", "NULL", "NULL")
        else:
            po = self.po(chrom, sv_left, sv_right, gene_name)
        return f"\t{self.sv(chrom, start, end)}\t{po}"
