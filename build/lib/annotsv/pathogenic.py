"""Pathogenic SV and pathogenic SNV/indel annotations: port of
``AnnotSV-pathogenicsv.tcl`` and ``AnnotSV-pathogenicsnvindel.tcl``.

* ``checkPathogenicFiles`` (+ the ClinVar, ClinGen, dbVar and OMIM morbid
  parsers)                    -> ``check_pathogenic_files``
* ``pathogenicSVannotation``   -> ``PathogenicSV``      (P_<type>_phen/hpo/source/coord)
* ``poPathogenicSVannotation`` -> ``PoPathogenicSV``    (po_P_<type>_phen/hpo/source/coord/percent)
* ``checkPathoSNVindelFile``   -> ``check_patho_snvindel_file``
* ``pathoSNVindelAnnotation``  -> ``PathoSNVindel``     (P_snvindel_nb, P_snvindel_phen)

``PathogenicAnnotations`` gathers the header and per-line text that
``OrganizeAnnotation`` writes right after ``RE_gene`` (Human only).

Tcl quirks kept on purpose (see tests/python/test_pathogenic_parity.py):

* ``checkPathogenicFiles`` loops over GRCh37, GRCh38 and CHM13 whatever the
  genome build, and then drops the ``P_<type>_*`` columns when the CHM13
  sorted BED is missing (``$pathogenicDir`` keeps its last loop value).
* The text of a line always gets "\\t" + the pathogenic text, even when no
  column was requested (an empty field with no header).
* poPathogenicSVannotation selects the svtypes with ``^po_P_`` for the empty
  text and the intersection, but with ``^P_`` when it builds the text of an
  annotated SV. Its ``po_P_*_hpo`` field is always empty (a misspelt
  variable in Tcl), and its coordinate regexp needs digits before ":", so
  pathogenic SVs on X/Y never count as partial overlaps.
* The empty P_snvindel text is "\\t" whatever the requested columns.
* checkPathoSNVindelFile fails (Tcl error, exit code 1) when the first ClinVar
  file it reads gives no SNV/indel.
"""

from __future__ import annotations

import glob
import gzip
import os
import re
import subprocess
from typing import Dict, List, Optional

from . import tcl
from .config import ConfigExit
from .gene_based import from_omim_gene_to_phenotype
from .general import check_bed, lines_from_file, write_text_in_file

SV_TYPES = ("Gain", "Loss", "Ins", "Inv")


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def _ann_dir(c) -> str:
    return f"{c['annotationsDir']}/Annotations_{c['organism']}"


def _patho_dir(c, build: str) -> str:
    return f"{_ann_dir(c)}/FtIncludedInSV/PathogenicSV/{build}"


def _snv_dir(c, build: str) -> str:
    return f"{_ann_dir(c)}/FtIncludedInSV/PathogenicSNVindel/{build}"


def _wanted(c) -> List[str]:
    cols = c["outputColHeader"]
    return tcl.split_list(cols) if isinstance(cols, str) else list(cols)


def _any_col(c, prefix_re: str) -> bool:
    """``lsearch -regexp $g_AnnotSV(outputColHeader) $prefix_re ne -1``."""
    return any(re.search(prefix_re, col) for col in _wanted(c))


def _gz_lines(path: str) -> List[str]:
    """``open "| gzip -cd $path"`` read with ``while {![eof $f]} {gets $f}``."""
    with gzip.open(path, "rb") as f:
        return f.read().decode("utf-8", errors="surrogateescape").split("\n")


def _read_lines(path: str) -> List[str]:
    with open(path, "r", encoding="utf-8", errors="surrogateescape", newline="") as f:
        return f.read().split("\n")


def _li(ls: List[str], i: int) -> str:
    return ls[i] if 0 <= i < len(ls) else ""


def _delete(path: str) -> None:
    if path and os.path.lexists(path):
        os.remove(path)


def _unix_sort(src: str, dst: str, caller: str) -> None:
    env = dict(os.environ, LC_ALL="C")
    with open(dst, "wb") as out:
        r = subprocess.run(["sort", "-k1,1", "-k2,2n", src], stdout=out, stderr=subprocess.PIPE, env=env)
    if r.returncode or r.stderr:
        raise ConfigExit(f"-- {caller} --\nsort -k1,1 -k2,2n {src} > {dst}\n"
                         f"{r.stderr.decode(errors='replace').rstrip(chr(10))}\nExit with error", 2)


def _intersect(c, a: str, out_file: str, full: bool, caller: str) -> None:
    """``bedtools intersect -sorted -a $a -b $fullAndSplitBedFile [-f 1.0] -wa -wb``,
    then without ``-sorted`` if it fails."""
    bt = c.get("bedtools", "bedtools")
    b = c["fullAndSplitBedFile"]
    f_opt = ["-f", "1.0"] if full else []

    def run(sorted_opt):
        with open(out_file, "w") as out:
            r = subprocess.run([bt, "intersect"] + sorted_opt + ["-a", a, "-b", b] + f_opt + ["-wa", "-wb"],
                               stdout=out, stderr=subprocess.PIPE)
        return r.returncode == 0 and not r.stderr, r.stderr.decode(errors="replace").rstrip("\n")

    ok, _ = run(["-sorted"])
    if ok:
        return
    ok, err = run([])
    if not ok:
        f_txt = " -f 1.0" if full else ""
        raise ConfigExit(f"-- {caller} --\n{bt} intersect -sorted -a {a} -b {b}{f_txt} -wa -wb > {out_file}\n"
                         f"{err}\nExit with error", 2)


def _out_tmp(c, repl: str) -> str:
    name = re.sub(r"(.formatted)?.bed$", lambda _m: repl, c["bedFile"], count=1, flags=re.I)
    return f"{c['outputDir']}/{os.path.basename(name)}"


# --------------------------------------------------------------------------
# ClinVar parsing (shared by the SV and SNV/indel files)
# --------------------------------------------------------------------------
_CLINVAR_SV_RE = re.compile(r"ALLELEID=([^;]+)?;.*CLNDISDB=([^;]+)?;.*CLNDN[INCL]*=([^;]+)?;.*"
                            r"CLNREVSTAT=([^;]+)?;.*CLNSIG=Pathogenic.*?CLNVC=([^;]+)?")
_CLINVAR_SNV_RE = re.compile(r"ALLELEID=([^;]+)?;.*CLNDISDB=([^;]+)?;.*CLNDN[INCL]*=([^;]+)?;.*"
                             r"CLNREVSTAT=([^;]+)?;.*CLNSIG=Pathogenic")
_HPO_RE = re.compile(r"Human_Phenotype_Ontology:HP:[0-9]+")


def _clinvar_records(path: str, regex, keep_size):
    """Yield (toWrite, CLNVC) for the selected ClinVar variants."""
    for L in _gz_lines(path):
        Ls = L.split("\t")
        if L.startswith("#"):
            continue
        m = regex.search(_li(Ls, 7))
        if not m:
            continue
        g = [x if x is not None else "" for x in m.groups()]
        alleleid, clndisdb, clndn, clnrevstat = g[0:4]
        clnvc = g[4] if len(g) > 4 else ""
        if not re.search("criteria_provided|_multiple_submitters|reviewed_by_expert_panel", clnrevstat):
            continue
        if "no_assertion_criteria_provided" in clnrevstat:
            continue
        ref, alt = _li(Ls, 3), _li(Ls, 4)
        sv_length = abs(len(ref) - len(alt)) + 1
        if not keep_size(sv_length):
            continue
        chrom, start = _li(Ls, 0), _li(Ls, 1)
        end = int(start) + sv_length
        coord = f"{chrom}:{start}-{end}"
        hpos = [h.replace("Human_Phenotype_Ontology:", "", 1) for h in (m2.group(0) for m2 in _HPO_RE.finditer(clndisdb))]
        clndn = re.sub(r"(\|)?not_provided", "", clndn, count=1)
        yield f"{chrom}\t{start}\t{end}\t{clndn}\t{tcl.format_list(hpos)}\tCLN:{alleleid}\t{coord}", clnvc


# --------------------------------------------------------------------------
# checkPathogenicFiles
# --------------------------------------------------------------------------
def _tmp_file(d: str, svtype: str, build: str) -> str:
    return f"{d}/pathogenic_{svtype}_SV_{build}.tmp.bed"


def _write(lines: List[str], path: str) -> None:
    if lines:
        write_text_in_file("\n".join(lines), path)


def _check_clinvar(c, build: str) -> None:
    d = _patho_dir(c, build)
    files = glob.glob(f"{d}/clinvar*.vcf.gz")
    if not files:
        return
    downloaded = files[-1]
    min_size = int(c["SVminSize"])
    to_write: Dict[str, List[str]] = {t: [] for t in SV_TYPES}
    kind = {"Deletion": "Loss", "Duplication": "Gain", "Insertion": "Ins", "Inversion": "Inv"}
    for text, clnvc in _clinvar_records(downloaded, _CLINVAR_SV_RE, lambda n: n >= min_size):
        if clnvc in kind:
            to_write[kind[clnvc]].append(text)
    for t in ("Loss", "Gain", "Ins", "Inv"):
        _write(to_write[t], _tmp_file(d, t, build))
    _delete(downloaded)


def _check_clingen(c, build: str) -> None:
    d = _patho_dir(c, build)
    f1 = f"{d}/ClinGen_gene_curation_list_{build}.tsv"
    f2 = f"{d}/ClinGen_region_curation_list_{build}.tsv"
    f1 = f1 if os.path.exists(f1) else ""
    f2 = f2 if os.path.exists(f2) else ""
    if f1 == "" and f2 == "":
        return
    loss: List[str] = []
    gain: List[str] = []
    for path, region in ((f1, False), (f2, True)):
        # Tcl: [open ""] fails when only one of the two files exists
        idx: Dict[str, int] = {}
        for L in _read_lines(path):
            Ls = L.split("\t")
            if re.search("(^#Gene Symbol)|(^#ISCA ID)", L):
                idx = {"id": 0}
                names = [("coord", "Genomic Location"), ("hi", "Haploinsufficiency Score"),
                         ("ts", "Triplosensitivity Score")]
                if region:
                    names.append(("isca", "ISCA Region Name"))
                for key, name in names:
                    if name not in Ls:
                        raise ConfigExit(f"\nBad header line syntax. \"{name}\" column not found - Exit with error", 2)
                    idx[key] = Ls.index(name)
                continue
            if L.startswith("#"):
                continue
            ID = _li(Ls, idx["id"])
            coord = re.sub("(chr)| ", "", _li(Ls, idx["coord"]))
            m = re.search(r"([0-9XYMT]+):([0-9]+)-([0-9]+)", coord)
            if not m:
                continue
            chrom, start, end = m.groups()
            hi, ts = _li(Ls, idx["hi"]), _li(Ls, idx["ts"])
            if hi == "3":
                phen = _li(Ls, idx["isca"]) if region else from_omim_gene_to_phenotype(c, ID)
                loss.append(f"{chrom}\t{start}\t{end}\t{phen}\t\tHI3:{ID}\t{coord}")
            if ts == "3":
                phen = _li(Ls, idx["isca"]) if region else from_omim_gene_to_phenotype(c, ID)
                gain.append(f"{chrom}\t{start}\t{end}\t{phen}\t\tTS3:{ID}\t{coord}")
    _write(loss, _tmp_file(d, "Loss", build))
    _write(gain, _tmp_file(d, "Gain", build))
    _delete(f1)
    _delete(f2)


def _check_dbvar(c, build: str) -> None:
    d = _patho_dir(c, build)
    files = []
    for kind in ("deletions", "duplications", "insertions"):
        p = f"{d}/{build}.nr_{kind}.pathogenic.tsv.gz"
        files.append(p if os.path.exists(p) else "")
    if not any(files):
        return
    to_write: Dict[str, List[str]] = {"Loss": [], "Gain": [], "Ins": []}
    svtype = None
    for path in files:
        if re.search("deletions", path):
            svtype = "Loss"
        elif re.search("duplications", path):
            svtype = "Gain"
        elif re.search("insertions", path):
            svtype = "Ins"
        if path == "":
            continue
        idx: Dict[str, int] = {}
        for L in _gz_lines(path):
            Ls = L.split("\t")
            if L.startswith("#chr"):
                idx = {"chr": 0}
                for key, name in (("start", "outermost_start"), ("stop", "outermost_stop"), ("source", "variant")):
                    if name not in Ls:
                        raise ConfigExit(f"{path}\nBad header line syntax. {name} column not found - Exit with error", 2)
                    idx[key] = Ls.index(name)
                continue
            if L[:1] == "#" or L == "":
                continue
            chrom, start, end = _li(Ls, idx["chr"]), _li(Ls, idx["start"]), _li(Ls, idx["stop"])
            to_write[svtype].append(f"{chrom}\t{start}\t{end}\t\t\tdbVar:{_li(Ls, idx['source'])}\t{chrom}:{start}-{end}")
    for t in ("Loss", "Gain", "Ins"):
        _write(to_write[t], _tmp_file(d, t, build))
    for p in files:
        _delete(p)


def _check_omim(c, build: str) -> None:
    d = _patho_dir(c, build)
    files = glob.glob(f"{d}/*_morbid.tsv.gz")
    if not files:
        return
    genes_bed = f"{_ann_dir(c)}/Genes/{build}/genes.{c['tx']}.sorted.bed"
    coords: Dict[str, tuple] = {}
    for L in lines_from_file(genes_bed if os.path.exists(genes_bed) else ""):
        Ls = L.split("\t")
        coords[_li(Ls, 4)] = (_li(Ls, 0), _li(Ls, 1), _li(Ls, 2))
    loss: List[str] = []
    for L in _gz_lines(files[0]):
        g = tcl.lindex(L, 0)
        if g == "":
            continue
        pheno = from_omim_gene_to_phenotype(c, g)
        if g in coords:
            chrom, start, end = coords[g]
            pheno = re.sub("{|}", "", pheno)
            loss.append(f"{chrom}\t{start}\t{end}\t{pheno}\t\tmorbid:{g}\t{chrom}:{start}-{end}")
    _write(loss, _tmp_file(d, "Loss", build))
    for p in files:
        _delete(p)


def check_pathogenic_files(config) -> None:
    """checkPathogenicFiles: build ``pathogenic_<type>_SV_<build>.sorted.bed``."""
    c = config
    d = ""
    for build in ("GRCh37", "GRCh38", "CHM13"):
        d = _patho_dir(c, build)
        for t in ("Loss", "Gain", "Ins", "Inv"):
            _delete(_tmp_file(d, t, build))
        _check_clinvar(c, build)
        _check_clingen(c, build)
        _check_dbvar(c, build)
        _check_omim(c, build)
        for t in ("Loss", "Gain", "Ins", "Inv"):
            sorted_bed = f"{d}/pathogenic_{t}_SV_{build}.sorted.bed"
            tmp = _tmp_file(d, t, build)
            if os.path.exists(tmp):
                if os.path.exists(sorted_bed):
                    # while {[gets $id_file1 line] >= 0} {puts $id_file2 $line}
                    with open(sorted_bed, "rb") as src:
                        lines = src.read().split(b"\n")
                    if lines[-1] == b"":
                        lines.pop()
                    with open(tmp, "ab") as dst:
                        dst.write(b"".join(line + b"\n" for line in lines))
                check_bed(tmp, c["outputDir"], d)
                _delete(tmp)
        for t in ("Loss", "Gain", "Ins", "Inv"):
            formatted = f"{d}/pathogenic_{t}_SV_{build}.tmp.formatted.bed"
            if not os.path.exists(formatted):
                continue
            _unix_sort(formatted, f"{d}/pathogenic_{t}_SV_{build}.sorted.bed", "checkPathogenicFiles")
            _delete(formatted)

    # Tcl: $pathogenicDir is the one of the last genome build (CHM13)
    build = "CHM13"
    for t in ("Gain", "Loss", "Ins", "Inv"):
        if not os.path.exists(f"{d}/pathogenic_{t}_SV_{build}.sorted.bed"):
            c["outputColHeader"] = [e for e in _wanted(c) if not re.search(f"^P_{t.lower()}_", e)]


# --------------------------------------------------------------------------
# pathogenicSVannotation / poPathogenicSVannotation
# --------------------------------------------------------------------------
class PathogenicSV:
    """pathogenicSVannotation: P_<type>_phen, _hpo, _source, _coord (pathogenic SV included in the SV)."""

    def __init__(self, config):
        self.config = config
        self.text: Optional[Dict[str, str]] = None
        self.empty = ""

    def _load(self) -> None:
        c = self.config
        build = c["genomeBuild"]
        d = _patho_dir(c, build)
        types = [t for t in SV_TYPES if _any_col(c, f"^P_{t.lower()}_")]
        self.empty = "\t".join("" for _ in types for _ in range(4))
        all_sv: List[str] = []
        phen: Dict[tuple, List[str]] = {}
        hpo: Dict[tuple, List[str]] = {}
        source: Dict[tuple, List[str]] = {}
        coord: Dict[tuple, List[str]] = {}
        for t in types:
            bed = f"{d}/pathogenic_{t}_SV_{build}.sorted.bed"
            bed = bed if os.path.exists(bed) else ""
            tmp = _out_tmp(c, f".intersect.pathogenic-{t}")
            _delete(tmp)
            _intersect(c, bed, tmp, True, f"pathogenicSVAnnotation, {t}")
            for L in _read_lines(tmp):
                if L == "":
                    continue
                Ls = L.split("\t")
                key = f"{_li(Ls, 7)},{_li(Ls, 8)},{_li(Ls, 9)}"
                all_sv.append(key)
                if _li(Ls, 3) != "":
                    phen.setdefault((key, t), []).append(_li(Ls, 3))
                if _li(Ls, 4) != "":
                    hpo.setdefault((key, t), []).extend(_li(Ls, 4).split(" "))
                source.setdefault((key, t), []).append(_li(Ls, 5))
                coord.setdefault((key, t), []).append(_li(Ls, 6))
            _delete(tmp)
        self.text = {}
        for key in sorted(set(all_sv)):
            parts: List[str] = []
            for t in types:
                k = (key, t)
                if k in coord:
                    parts.append(";".join(sorted(set(phen[k]))) if k in phen else "")
                    parts.append(";".join(sorted(set(hpo[k]))) if k in hpo else "")
                    parts.append(";".join(sorted(set(source[k]))))
                    parts.append(";".join(sorted(set(coord[k]))))
                else:
                    parts += ["", "", "", ""]
            self.text[key] = "\t".join(parts)

    def __call__(self, chrom: str, start: str, end: str) -> str:
        if self.text is None:
            self._load()
        return self.text.get(f"{chrom},{start},{end}", self.empty)


_PO_COORD_RE = re.compile(r"[0-9]+:([0-9]+)-([0-9]+)")


class PoPathogenicSV:
    """poPathogenicSVannotation: po_P_<gain|loss>_phen, _hpo, _source, _coord, _percent."""

    def __init__(self, config):
        self.config = config
        self.text: Optional[Dict[str, str]] = None
        self.empty = ""

    def _load(self) -> None:
        c = self.config
        build = c["genomeBuild"]
        d = _patho_dir(c, build)
        min_size = int(c["SVminSize"])
        types = [t for t in ("Gain", "Loss") if _any_col(c, f"^po_P_{t.lower()}_")]
        self.empty = "\t".join("" for _ in types for _ in range(5))
        all_sv: List[str] = []
        phen: Dict[tuple, List[str]] = {}
        source: Dict[tuple, List[str]] = {}
        coord: Dict[tuple, List[str]] = {}
        percent: Dict[tuple, List[str]] = {}
        for t in types:
            bed = f"{d}/pathogenic_{t}_SV_{build}.sorted.bed"
            bed = bed if os.path.exists(bed) else ""
            tmp = _out_tmp(c, f".intersect.po_pathogenic-{t}")
            _delete(tmp)
            _intersect(c, bed, tmp, False, f"poPathogenicSVAnnotation, {t}")
            for L in _read_lines(tmp):
                if L == "":
                    continue
                Ls = L.split("\t")
                m = _PO_COORD_RE.search(_li(Ls, 6))
                if not m:
                    continue
                start_p, end_p = int(m.group(1)), int(m.group(2))
                if end_p - start_p < min_size:
                    continue
                sv_start, sv_end = int(_li(Ls, 8)), int(_li(Ls, 9))
                po_end = end_p if sv_end > end_p else sv_end
                po_start = sv_start if sv_start > start_p else start_p
                length = end_p - start_p
                if length == 0:
                    length = 1
                pct = "%.2f" % ((po_end - po_start) * 100.0 / length)
                if float(pct) == 100:
                    continue
                key = f"{_li(Ls, 7)},{_li(Ls, 8)},{_li(Ls, 9)}"
                all_sv.append(key)
                k = (key, t)
                if _li(Ls, 3) != "":
                    phen.setdefault(k, []).append(_li(Ls, 3))
                source.setdefault(k, []).append(_li(Ls, 5))
                coord.setdefault(k, []).append(_li(Ls, 6))
                percent.setdefault(k, []).append(pct)
            _delete(tmp)
        self.text = {}
        # Tcl quirk: "^P_" (not "^po_P_") selects the svtypes of the annotated SVs
        text_types = [t for t in ("Gain", "Loss") if _any_col(c, f"^P_{t.lower()}_")]
        for key in sorted(set(all_sv)):
            parts: List[str] = []
            for t in text_types:
                k = (key, t)
                if k in coord:
                    parts.append(";".join(sorted(set(phen[k]))) if k in phen else "")
                    parts.append("")  # Tcl quirk: misspelt hpo variable
                    parts.append(";".join(source[k]))
                    parts.append(";".join(coord[k]))
                    parts.append(";".join(percent[k]))
                else:
                    parts += ["", "", "", "", ""]
            self.text[key] = "\t".join(parts)

    def __call__(self, chrom: str, start: str, end: str) -> str:
        if self.text is None:
            self._load()
        return self.text.get(f"{chrom},{start},{end}", self.empty)


# --------------------------------------------------------------------------
# checkPathoSNVindelFile / pathoSNVindelAnnotation
# --------------------------------------------------------------------------
def check_patho_snvindel_file(config) -> None:
    """checkPathoSNVindelFile: build ``pathogenic_SNVindel_<build>.sorted.bed``."""
    c = config
    min_size = int(c["SVminSize"])
    formatted: Optional[str] = None  # Tcl variable kept across the genome builds
    for build in ("GRCh37", "GRCh38"):
        d = _snv_dir(c, build)
        files = glob.glob(f"{d}/clinvar*.vcf.gz")
        if not files:
            continue
        downloaded = files[-1]
        bed = f"{d}/pathogenic_SNVindel_{build}.bed"
        _delete(bed)
        to_write = [text for text, _ in _clinvar_records(downloaded, _CLINVAR_SNV_RE, lambda n: n < min_size)]
        _write(to_write, bed)
        _delete(downloaded)
        if os.path.exists(bed):
            formatted = re.sub(r".bed$", ".formatted.bed", bed, count=1)
            check_bed(bed, c["outputDir"], d)
            _delete(bed)
        if formatted is None:
            # Tcl error (tclsh exits with 1)
            raise ConfigExit("can't read \"pathoSNVindelBEDfile_Formatted\": no such variable", 1)
        if os.path.exists(formatted):
            sorted_bed = re.sub(r".bed$", ".sorted.bed", bed, count=1)
            _unix_sort(formatted, sorted_bed, "checkPathoSNVindelFile")
            _delete(formatted)


class PathoSNVindel:
    """pathoSNVindelAnnotation: "P_snvindel_nb\\tP_snvindel_phen"."""

    def __init__(self, config):
        self.config = config
        self.text: Optional[Dict[str, str]] = None
        self.empty = "\t"

    def _load(self) -> None:
        c = self.config
        build = c["genomeBuild"]
        bed = f"{_snv_dir(c, build)}/pathogenic_SNVindel_{build}.sorted.bed"
        self.text = {}
        if not os.path.exists(bed):
            return
        tmp = _out_tmp(c, ".intersect.pathogenic-snvindel")
        _delete(tmp)
        _intersect(c, bed, tmp, True, "pathoSNVindelAnnotation")
        nb: Dict[str, int] = {}
        phen: Dict[str, List[str]] = {}
        for L in _read_lines(tmp):
            if L == "":
                continue
            Ls = L.split("\t")
            key = f"{_li(Ls, 7)},{_li(Ls, 8)},{_li(Ls, 9)}"
            lst = phen.setdefault(key, [])
            lst += [p for p in _li(Ls, 3).split("|") if p != ""]
            nb[key] = nb.get(key, 0) + 1
        _delete(tmp)
        for key in nb:
            self.text[key] = f"{nb[key]}\t" + ";".join(sorted(set(phen[key])))

    def __call__(self, chrom: str, start: str, end: str) -> str:
        if self.text is None:
            self._load()
        return self.text.get(f"{chrom},{start},{end}", self.empty)


# --------------------------------------------------------------------------
# OrganizeAnnotation part
# --------------------------------------------------------------------------
class PathogenicAnnotations:
    """The pathogenic columns that OrganizeAnnotation writes after ``RE_gene`` (Human only)."""

    def __init__(self, config):
        self.config = config
        self.sv = PathogenicSV(config)
        self.po = PoPathogenicSV(config)
        self.snv = PathoSNVindel(config)

    @property
    def human(self) -> bool:
        return self.config["organism"] == "Human"

    def header(self) -> str:
        c = self.config
        if not self.human:
            return ""
        h = ""
        for t in ("gain", "loss", "ins", "inv"):
            if _any_col(c, f"^P_{t}_"):
                h += f"\tP_{t}_phen\tP_{t}_hpo\tP_{t}_source\tP_{t}_coord"
        for t in ("gain", "loss"):
            if _any_col(c, f"^po_P_{t}_"):
                h += f"\tpo_P_{t}_phen\tpo_P_{t}_hpo\tpo_P_{t}_source\tpo_P_{t}_coord\tpo_P_{t}_percent"
        if _any_col(c, "^P_snvindel_"):
            h += "\tP_snvindel_nb\tP_snvindel_phen"
        return h

    def text(self, chrom: str, start: str, end: str, sv_left: str, sv_right: str, mode: str) -> str:
        """``start``/``end``: the intersection on split lines, the SV on full lines."""
        if not self.human:
            return ""
        if mode == "split":
            po = self.po("NULL", "NULL", "NULL")
        else:
            po = self.po(chrom, sv_left, sv_right)
        return f"\t{self.sv(chrom, start, end)}\t{po}\t{self.snv(chrom, start, end)}"
