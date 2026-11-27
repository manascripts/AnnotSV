"""Regulatory elements: port of ``AnnotSV-regulatoryelements.tcl``.

Two parts:

* the ``check*`` procs that build the formatted ``*.sorted.bed`` files of
  ``FtIncludedInSV/RegulatoryElements`` from the downloaded sources
  (promoters, EnhancerAtlas, GeneHancer, miRTargetLink, ABC and MPRA), once;
* ``regulatory_elements_annotation`` (``regulatoryElementsAnnotation``), which
  intersects the SVs with these files and builds the ``RE_gene`` text of each
  SV (``g_re``) and the ClinGen HI/TS of the regulated genes (``g_HI``,
  ``g_TS``).

The Tcl quirks are kept (see tests/reference/PARITY.md): the ``glob -no
complain`` of checkGHfiles, ``lindex`` on whole lines, stale values after a
badly formatted MPRA line, the blank lines that empty chromosomes leave in the
formatted files, the trailing empty line that ``while {![eof $f]}`` loops read,
``PG=A`` written for any PhenoGenius specificity when ``-REselect1 0``, etc.
"""

from __future__ import annotations

import glob
import gzip
import os
import re
import subprocess
import sys
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional

from . import gene_based, tcl
from .config import ConfigExit
from .general import write_text_in_file
from .genes import genes_dir

CHROMS = "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M MT".split()


# --------------------------------------------------------------------------
# Tcl file reading helpers (channels use the "auto" end of line translation)
# --------------------------------------------------------------------------
def _read(path: str) -> str:
    if path.endswith(".gz"):
        with gzip.open(path, "rt", errors="surrogateescape") as f:
            return f.read()
    with open(path, "r", errors="surrogateescape") as f:
        return f.read()


def _lines_from_file(path: str) -> List[str]:
    """LinesFromFile: ``split [read -nonewline] "\\n"``."""
    text = _read(path)
    return (text[:-1] if text.endswith("\n") else text).split("\n")


def _eof_loop_lines(text: str) -> List[str]:
    """Lines read by ``while {![eof $f]} {set L [gets $f] ...}``: an empty line
    comes last when the text ends with a newline (or is empty)."""
    return text.split("\n")


def _gets_lines(text: str) -> List[str]:
    """Lines read by ``while {[gets $f L] >= 0}``."""
    if text == "":
        return []
    return (text[:-1] if text.endswith("\n") else text).split("\n")


def _tsplit(line: str) -> List[str]:
    """``split $L "\\t"`` (an empty string gives an empty list)."""
    return line.split("\t") if line != "" else []


def _li(ls: List[str], i: Optional[int]) -> str:
    """``lindex $Ls $i`` on a list (out of range, or ``-1`` from lsearch: "")."""
    if i is None:
        raise ConfigExit('can\'t read "i_...": no such variable', 1)
    return ls[i] if 0 <= i < len(ls) else ""


def _lsearch(ls: List[str], value: str) -> int:
    return ls.index(value) if value in ls else -1


def _now() -> str:
    return time.strftime("%B %d %Y - %H:%M")


def _sort_lines(lines: List[str]) -> List[str]:
    """``lsort -command AscendingSortOnElement1 [lsort -command AscendingSortOnElement2 $lines]``."""
    return tcl.sort_on_element(tcl.sort_on_element(lines, 2), 1)


def _rm(path: str) -> None:
    if os.path.lexists(path):
        if os.path.isdir(path) and not os.path.islink(path):
            import shutil
            shutil.rmtree(path)
        else:
            os.remove(path)


def _sort_bed(src: str, dst: str, caller: str, redirect: str = ">>") -> None:
    """The sort.tmp.bash script: ``LC_ALL=C sort -k1,1 -k2,2n src >> dst``.

    As with Tcl ``exec``, any stderr output is a failure (exit 2). The
    redirection creates ``dst`` even when sort fails.
    """
    env = dict(os.environ, LC_ALL="C")
    with open(dst, "ab" if redirect == ">>" else "wb") as out:
        r = subprocess.run(["sort", "-k1,1", "-k2,2n", src], stdout=out, stderr=subprocess.PIPE, env=env)
    if r.returncode or r.stderr:
        msg = r.stderr.decode(errors="replace").rstrip("\n")
        if not msg:
            msg = "child process exited abnormally"
        raise ConfigExit(f"-- {caller} --\nsort -k1,1 -k2,2n {src} {redirect} {dst}\n{msg}\nExit with error", 2)
    _rm(src)


def _reg_dir(config) -> str:
    c = config
    return f"{c['annotationsDir']}/Annotations_{c['organism']}/FtIncludedInSV/RegulatoryElements"


# --------------------------------------------------------------------------
# Gene name lists (isARefSeqGeneName, isAnENSEMBLgeneName, isMorbid, isCandidate)
# --------------------------------------------------------------------------
class GeneLists:
    """The lazily loaded gene sets of ``g_exist`` (``g_morbid`` and ``g_candidate`` are in gene_based)."""

    def __init__(self, config):
        self.config = config
        self._refseq = None
        self._ensembl = None

    def is_refseq(self, gene: str) -> bool:
        if self._refseq is None:
            path = f"{genes_dir(self.config)}/genes.RefSeq.sorted.bed"
            self._refseq = {_li(_tsplit(L), 4) for L in _lines_from_file(path)}
        return gene in self._refseq

    def is_ensembl(self, gene: str) -> bool:
        if self._ensembl is None:
            path = f"{genes_dir(self.config)}/genes.ENSEMBL.sorted.bed"
            self._ensembl = set()
            if os.path.exists(path):
                self._ensembl = {_li(_tsplit(L), 4) for L in _lines_from_file(path)}
        return gene in self._ensembl

    def is_morbid(self, gene: str) -> bool:
        """isMorbid (AnnotSV-omim.tcl, ported in gene_based)."""
        return gene_based.is_morbid(self.config, gene)

    def is_candidate(self, gene: str) -> bool:
        """isCandidate (AnnotSV-candidategenes.tcl, ported in gene_based)."""
        return gene_based.is_candidate(self.config, gene)


def _gene_lists(config) -> GeneLists:
    lists = getattr(config, "_gene_lists", None)
    if lists is None:
        lists = GeneLists(config)
        try:
            setattr(config, "_gene_lists", lists)
        except AttributeError:
            pass
    return lists


# --------------------------------------------------------------------------
# checkPromoterFile
# --------------------------------------------------------------------------
def check_promoter_file(config) -> None:
    c = config
    c["promAnn"] = 0
    reg_dir = f"{_reg_dir(c)}/{c['genomeBuild']}"
    os.makedirs(reg_dir, exist_ok=True)
    genes_file = f"{genes_dir(c)}/genes.{c['tx']}.sorted.bed"
    formatted = f"{reg_dir}/promoter_{c['promoterSize']}bp_{c['tx']}_{c['genomeBuild']}.sorted.bed"

    if not os.path.exists(formatted):
        print(f"\t...{c['tx']} promoters configuration ({_now()})")
        print(f"\t   (--promoter-size = {c['promoterSize']})\n")
        size = int(c["promoterSize"])
        genes_prom: Dict[str, List[str]] = {}
        l_prom: Dict[str, List[str]] = {}
        prom_left = prom_right = None
        for L in _lines_from_file(genes_file):
            ls = _tsplit(L)
            chrom = _li(ls, 0)
            strand = tcl.lindex(L, 3)
            if strand == "+":
                prom_right = int(_li(ls, 1))
                prom_left = max(prom_right - size, 0)
            elif strand == "-":
                prom_left = int(_li(ls, 2))
                prom_right = prom_left + size
            else:
                print(f'strand eq "{strand}" -- {L}')
                if prom_left is None:
                    raise ConfigExit("can't read \"promLeft\": no such variable", 1)
            genes_prom.setdefault(f"{chrom}\t{prom_left}\t{prom_right}", []).append(_li(ls, 4))
            l_prom.setdefault(chrom, []).append(f"{prom_left}\t{prom_right}")
        all_prom = []
        for chrom in CHROMS:
            if chrom not in l_prom:
                continue
            proms = tcl.sort_on_element(tcl.sort_on_element(sorted(set(l_prom[chrom])), 1), 0)
            for prom in proms:
                genes = ";".join(sorted(set(genes_prom[f"{chrom}\t{prom}"])))
                all_prom.append(f"{chrom}\t{prom}\t{c['tx']}_promoter\t{genes}")
        write_text_in_file("\n".join(all_prom), f"{formatted}.tmp")
        _sort_bed(f"{formatted}.tmp", formatted, "checkPromoterFile", ">")
    c["promAnn"] = 1


# --------------------------------------------------------------------------
# EnhancerAtlas and ABC (same layout: one coordinate -> genes file per tx)
# --------------------------------------------------------------------------
def _write_coord_genes(coords: Dict[str, List[str]], genes: Dict[str, List[str]], label: str, out: str) -> None:
    for chrom in CHROMS:
        if chrom not in coords:
            continue
        lines = []
        for coord in _sort_lines(sorted(set(coords[chrom]))):
            lines.append(f"{coord}\t{label}\t{';'.join(sorted(set(genes[coord])))}")
        write_text_in_file("\n".join(lines), out)


_EA_RE = re.compile(r"chr([0-9XYMT]+):([0-9]+)-([0-9]+)_[^$]+\$([^$]+)")


def check_ea_files(config) -> None:
    c = config
    c["EAann"] = 0
    reg_dir = _reg_dir(c)
    build = c["genomeBuild"]
    refseq_f = f"{reg_dir}/{build}/EA_RefSeq_{build}.sorted.bed"
    ensembl_f = f"{reg_dir}/{build}/EA_ENSEMBL_{build}.sorted.bed"
    necessary = f"{reg_dir}/{build}/EA_{c['tx']}_{build}.sorted.bed"
    if os.path.exists(necessary):
        c["EAann"] = 1
    if os.path.exists(refseq_f) and os.path.exists(ensembl_f):
        return
    ea_files = glob.glob(f"{reg_dir}/GRCh37/*_EP.txt")
    if not ea_files:
        return

    print(f"\t...EnhancerAtlas configuration ({_now()})")
    print(f"\t\t...creation of {os.path.basename(refseq_f)} and {os.path.basename(ensembl_f)}")
    print("\t\t   (done only once)")
    print("\t\t   For GRCh38 use, please lift over these GRCh37 files to GRCh38")
    for p in (refseq_f, ensembl_f, f"{refseq_f}.tmp", f"{ensembl_f}.tmp"):
        _rm(p)

    gl = _gene_lists(c)
    coord_r: Dict[str, List[str]] = {}
    genes_r: Dict[str, List[str]] = {}
    coord_e: Dict[str, List[str]] = {}
    genes_e: Dict[str, List[str]] = {}
    for ea_file in ea_files:
        for L in _lines_from_file(ea_file):
            m = _EA_RE.search(L)
            if not m:
                continue
            chrom, start, end, gene = m.groups()
            coord = f"{chrom}\t{start}\t{end}"
            if gl.is_refseq(gene):
                coord_r.setdefault(chrom, []).append(coord)
                genes_r.setdefault(coord, []).append(gene)
            if gl.is_ensembl(gene):
                coord_e.setdefault(chrom, []).append(coord)
                genes_e.setdefault(coord, []).append(gene)
    _write_coord_genes(coord_r, genes_r, "EA_enhancer", f"{refseq_f}.tmp")
    _write_coord_genes(coord_e, genes_e, "EA_enhancer", f"{ensembl_f}.tmp")
    _sort_bed(f"{refseq_f}.tmp", refseq_f, "checkEAfiles")
    _sort_bed(f"{ensembl_f}.tmp", ensembl_f, "checkEAfiles")

    if os.path.exists(refseq_f) and os.path.exists(ensembl_f):
        for p in ea_files:
            _rm(p)
    if os.path.exists(necessary):
        c["EAann"] = 1


def check_abc_files(config) -> None:
    c = config
    c["ABCann"] = 0
    reg_dir = _reg_dir(c)
    build = c["genomeBuild"]
    necessary = f"{reg_dir}/{build}/ABC_{c['tx']}_{build}.sorted.bed"
    if os.path.exists(necessary):
        c["ABCann"] = 1
    refseq_f = f"{reg_dir}/GRCh37/ABC_RefSeq_GRCh37.sorted.bed"
    ensembl_f = f"{reg_dir}/GRCh37/ABC_ENSEMBL_GRCh37.sorted.bed"
    if os.path.exists(refseq_f) and os.path.exists(ensembl_f):
        return
    downloaded = f"{reg_dir}/GRCh37/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz"
    if not os.path.exists(downloaded):
        return

    print(f"\t...ABC configuration ({_now()})")
    print(f"\t\t...creation of {os.path.basename(refseq_f)} and {os.path.basename(ensembl_f)}")
    print("\t\t   (done only once)")
    print("\t\t   For GRCh38 use, please lift over these GRCh37 files to GRCh38")
    for p in (refseq_f, ensembl_f, f"{refseq_f}.tmp", f"{ensembl_f}.tmp"):
        _rm(p)

    gl = _gene_lists(c)
    coord_r: Dict[str, List[str]] = {}
    genes_r: Dict[str, List[str]] = {}
    coord_e: Dict[str, List[str]] = {}
    genes_e: Dict[str, List[str]] = {}
    i_start = i_end = i_gene = None
    for L in _eof_loop_lines(_read(downloaded)):
        ls = _tsplit(L)
        if L.startswith("chr\t"):
            i_start, i_end, i_gene = _lsearch(ls, "start"), _lsearch(ls, "end"), _lsearch(ls, "TargetGene")
            continue
        chrom = _li(ls, 0).replace("chr", "", 1)
        start, end, gene = _li(ls, i_start), _li(ls, i_end), _li(ls, i_gene)
        coord = f"{chrom}\t{start}\t{end}"
        if gl.is_refseq(gene):
            coord_r.setdefault(chrom, []).append(coord)
            genes_r.setdefault(coord, []).append(gene)
        if gl.is_ensembl(gene):
            coord_e.setdefault(chrom, []).append(coord)
            genes_e.setdefault(coord, []).append(gene)
    _write_coord_genes(coord_r, genes_r, "ABC_enhancer", f"{refseq_f}.tmp")
    _write_coord_genes(coord_e, genes_e, "ABC_enhancer", f"{ensembl_f}.tmp")
    _sort_bed(f"{refseq_f}.tmp", refseq_f, "checkABCfiles")
    _sort_bed(f"{ensembl_f}.tmp", ensembl_f, "checkABCfiles")

    if os.path.exists(refseq_f) and os.path.exists(ensembl_f):
        _rm(downloaded)
    if os.path.exists(necessary):
        c["ABCann"] = 1


# --------------------------------------------------------------------------
# checkGHfiles
# --------------------------------------------------------------------------
def _glob_last(patterns: List[str]) -> str:
    found: List[str] = []
    for p in patterns:
        found += glob.glob(p)
    return found[-1] if found else "no_file_exists"


def check_gh_files(config) -> None:
    c = config
    c["GHann"] = 0
    reg_dir = _reg_dir(c)
    build = c["genomeBuild"]
    if (os.path.exists(f"{reg_dir}/{build}/GH_RefSeq_{build}.sorted.bed")
            and os.path.exists(f"{reg_dir}/{build}/GH_ENSEMBL_{build}.sorted.bed")):
        c["GHann"] = 1
        return

    # Tcl: [glob -no complain "..."]: "-no" stands for -nocomplain and "complain" is one more pattern
    elements_f = _glob_last(["complain", f"{reg_dir}/GeneHancer*elements*.txt"])
    associations_f = _glob_last([f"{reg_dir}/GeneHancer*gene_association*_scores*.txt"])
    hg19_f = _glob_last([f"{reg_dir}/GeneHancer*hg19*.txt"])
    for p in (elements_f, associations_f, hg19_f):
        if not os.path.exists(p):
            if c["organism"] == "Human":
                print("\nWARNING: No GeneHancer annotations available.")
                print("(Please, see in the README file how to add these annotations. "
                      "Users need to contact the GeneCards team.)\n")
            return

    print(f"\t...GeneHancer configuration ({_now()})")
    out = {
        ("RefSeq", "37"): f"{reg_dir}/GRCh37/GH_RefSeq_GRCh37.sorted.bed",
        ("RefSeq", "38"): f"{reg_dir}/GRCh38/GH_RefSeq_GRCh38.sorted.bed",
        ("ENSEMBL", "37"): f"{reg_dir}/GRCh37/GH_ENSEMBL_GRCh37.sorted.bed",
        ("ENSEMBL", "38"): f"{reg_dir}/GRCh38/GH_ENSEMBL_GRCh38.sorted.bed",
    }
    print("\t\t...creation of the 4 GH_*_GRCh3*.sorted.bed files")
    print("\t\t   (done only once)")
    for p in out.values():
        _rm(p)

    gl = _gene_lists(c)
    l_ghid: Dict[str, List[str]] = {}
    coord38: Dict[str, str] = {}
    coord37: Dict[str, str] = {}
    gh_type: Dict[str, str] = {}
    gh_gene: Dict[tuple, List[str]] = {}
    i_start = i_end = i_id = i_type = None
    for L in _lines_from_file(elements_f):
        ls = _tsplit(L)
        if L.startswith("chr"):
            i_start, i_end = _lsearch(ls, "element_start"), _lsearch(ls, "element_end")
            i_id, i_type = _lsearch(ls, "GHid"), _lsearch(ls, "regulatory_element_type")
            continue
        chrom = _li(ls, 0).replace("chr", "", 1)
        ghid = _li(ls, i_id)
        l_ghid.setdefault(chrom, []).append(ghid)
        coord38[ghid] = f"{chrom}\t{_li(ls, i_start)}\t{_li(ls, i_end)}"
        t = _li(ls, i_type).replace("/", "_").replace("promoter_enhancer", "enhancer", 1)
        gh_type[ghid] = f"GH_{t.lower()}"
        gh_gene[("RefSeq", ghid)] = []
        gh_gene[("ENSEMBL", ghid)] = []

    i_id = i_gene = None
    for L in _lines_from_file(associations_f):
        ls = _tsplit(L)
        if L.startswith("GHid"):
            i_id, i_gene = _lsearch(ls, "GHid"), _lsearch(ls, "symbol")
            continue
        ghid, gene = _li(ls, i_id), _li(ls, i_gene)
        if gl.is_refseq(gene):
            gh_gene.setdefault(("RefSeq", ghid), []).append(gene)
        if gl.is_ensembl(gene):
            gh_gene.setdefault(("ENSEMBL", ghid), []).append(gene)

    i_start = i_end = i_id = None
    for L in _lines_from_file(hg19_f):
        ls = _tsplit(L)
        if L.startswith("chromosome"):
            i_start, i_end, i_id = _lsearch(ls, "start"), _lsearch(ls, "end"), _lsearch(ls, "GHid")
            continue
        ghid = _li(ls, i_id)
        chrom = _li(ls, 0).replace("chr", "", 1)
        coord37[ghid] = f"{chrom}\t{_li(ls, i_start)}\t{_li(ls, i_end)}"

    for p in out.values():
        _rm(f"{p}.tmp")
    for chrom in CHROMS:
        if chrom not in l_ghid:
            continue
        lines: Dict[tuple, List[str]] = {k: [] for k in out}
        for ghid in l_ghid[chrom]:
            for tx in ("RefSeq", "ENSEMBL"):
                genes = gh_gene[(tx, ghid)]
                if not genes:
                    continue
                for b, coords in (("37", coord37), ("38", coord38)):
                    if ghid in coords:  # catch {...}: no line without coordinates
                        lines[(tx, b)].append(f"{coords[ghid]}\t{gh_type[ghid]}\t{';'.join(genes)}")
        for k in (("RefSeq", "37"), ("RefSeq", "38"), ("ENSEMBL", "37"), ("ENSEMBL", "38")):
            write_text_in_file("\n".join(_sort_lines(lines[k])), f"{out[k]}.tmp")

    for k in (("RefSeq", "37"), ("RefSeq", "38"), ("ENSEMBL", "37"), ("ENSEMBL", "38")):
        _sort_bed(f"{out[k]}.tmp", out[k], "checkGHfiles")

    if all(os.path.exists(p) for p in out.values()):
        for p in (elements_f, associations_f, hg19_f, f"{reg_dir}/GeneHancer_tissues.txt", f"{reg_dir}/ReadMe.txt"):
            _rm(p)
    c["GHann"] = 1


# --------------------------------------------------------------------------
# checkMiRTargetLinkFiles
# --------------------------------------------------------------------------
def check_mirtargetlink_files(config) -> None:
    c = config
    c["miRNAann"] = 0
    reg_dir = _reg_dir(c)
    build = c["genomeBuild"]
    refseq_f = f"{reg_dir}/{build}/miRTargetLink_RefSeq_{build}.sorted.bed"
    ensembl_f = f"{reg_dir}/{build}/miRTargetLink_ENSEMBL_{build}.sorted.bed"
    if c["organism"] == "Human":
        if os.path.exists(refseq_f) and os.path.exists(ensembl_f):
            c["miRNAann"] = 1
            return
    elif c["organism"] == "Mouse":
        if os.path.exists(refseq_f):
            c["miRNAann"] = 1
            return
    if build not in ("GRCh38", "mm10"):
        return
    if build == "GRCh38":
        downloaded = glob.glob(f"{reg_dir}/*_miRNA-gene_pairs_hsa_miRBase_*_GRCh38_location_augmented.tsv")
    else:
        downloaded = glob.glob(f"{reg_dir}/*_miRNA-gene_pairs_mmu_miRBase_*_GRCm38_location_augmented.tsv")
    if not downloaded:
        return

    print(f"\t...miRTargetLink configuration ({_now()})")
    print(f"\t\t...creation of the miRNA_*_{build}.sorted.bed files")
    print("\t\t   (done only once)")
    _rm(refseq_f)
    _rm(ensembl_f)

    l_genes: Dict[str, List[str]] = {}
    l_coord: Dict[str, List[str]] = {}
    i_chrom = i_start = i_end = i_gene = None
    for mirna_file in downloaded:
        for L in _lines_from_file(mirna_file):
            ls = _tsplit(L)
            if L.startswith("seqnames"):
                i_chrom, i_start = _lsearch(ls, "seqnames"), _lsearch(ls, "start")
                i_end, i_gene = _lsearch(ls, "end"), _lsearch(ls, "target_gene")
                continue
            chrom = _li(ls, i_chrom).replace("chr", "", 1)
            coord = f"{chrom}\t{_li(ls, i_start)}\t{_li(ls, i_end)}"
            l_genes.setdefault(coord, []).append(_li(ls, i_gene))
            l_coord.setdefault(chrom, []).append(coord)

    gl = _gene_lists(c)
    _rm(f"{refseq_f}.tmp")
    _rm(f"{ensembl_f}.tmp")
    for chrom in CHROMS:
        if chrom not in l_coord:
            continue
        lines_r: List[str] = []
        lines_e: List[str] = []
        for coord in sorted(set(l_coord[chrom])):
            g_r = [g for g in l_genes[coord] if gl.is_refseq(g)]
            g_e = [g for g in l_genes[coord] if gl.is_ensembl(g)]
            if g_r:
                lines_r.append(f"{coord}\tmTL_miRNA\t{';'.join(sorted(set(g_r)))}")
            if g_e:
                lines_e.append(f"{coord}\tmTL_miRNA\t{';'.join(sorted(set(g_e)))}")
        write_text_in_file("\n".join(_sort_lines(lines_r)), f"{refseq_f}.tmp")
        if build == "GRCh38":
            write_text_in_file("\n".join(_sort_lines(lines_e)), f"{ensembl_f}.tmp")

    for p in (refseq_f, ensembl_f):
        if os.path.exists(f"{p}.tmp"):
            _sort_bed(f"{p}.tmp", p, "checkMiRTargetLinkFiles")
    for p in downloaded:
        _rm(p)
    c["miRNAann"] = 1


# --------------------------------------------------------------------------
# checkMPRAfiles
# --------------------------------------------------------------------------
_MPRA_COORD = re.compile(r"chr([0-9XYMT]+):([0-9]+)-([0-9]+)")


def check_mpra_files(config) -> None:
    c = config
    c["MPRAann"] = 0
    reg_dir = _reg_dir(c)
    build = c["genomeBuild"]
    necessary = f"{reg_dir}/{build}/MPRA_{c['tx']}_{build}.sorted.bed"
    if os.path.exists(necessary):
        c["MPRAann"] = 1
        return
    out = {(tx, b): f"{reg_dir}/GRCh{b}/MPRA_{tx}_GRCh{b}.sorted.bed"
           for tx in ("RefSeq", "ENSEMBL") for b in ("37", "38")}
    if all(os.path.exists(p) for p in out.values()):
        return
    prom_f = f"{reg_dir}/MPRA_promoters.tsv"
    enh_f = f"{reg_dir}/MPRA_enhancers.tsv"
    if not os.path.exists(prom_f) or not os.path.exists(enh_f):
        return

    print(f"\t...MPRA configuration ({_now()})")
    print(f"\t\t...creation of {os.path.basename(out[('RefSeq', '37')])}, {os.path.basename(out[('ENSEMBL', '37')])}, "
          f"{os.path.basename(out[('RefSeq', '38')])} and {os.path.basename(out[('ENSEMBL', '38')])}")
    print("\t\t   (done only once)")
    for p in out.values():
        _rm(p)
    for p in out.values():
        _rm(f"{p}.tmp")

    gl = _gene_lists(c)
    mpra_type: Dict[str, str] = {}
    coords: Dict[tuple, Dict[str, List[str]]] = {k: {} for k in out}
    genes: Dict[tuple, Dict[str, List[str]]] = {k: {} for k in out}
    i37 = i38 = None
    t = None
    m37 = m38 = None  # the last coordinates that matched (Tcl keeps them after a bad line)
    for L in _lines_from_file(prom_f) + _lines_from_file(enh_f):
        ls = _tsplit(L)
        if L.startswith("Name\t"):
            i37 = _lsearch(ls, "Genomic coordinates (GRCh37)")
            i38 = _lsearch(ls, "Genomic coordinates (GRCh38)")
            t = "MPRA_promoter" if "Transcript" in L else "MPRA_enhancer"
            continue
        c37 = _li(ls, i37).replace(",", "")
        c38 = _li(ls, i38).replace(",", "")
        for text, which in ((c37, 37), (c38, 38)):
            m = _MPRA_COORD.search(text)
            if m:
                if which == 37:
                    m37 = m.groups()
                else:
                    m38 = m.groups()
            else:
                print(f"Bad format in {prom_f}")
                print(L)
                print("Check this format.")
        if m37 is None or m38 is None or t is None:
            raise ConfigExit("can't read \"chrom37\": no such variable", 1)
        gene = re.sub(r" \(.+|\+.+", "", _li(ls, 0), count=1)
        k37, k38 = "\t".join(m37), "\t".join(m38)
        mpra_type[k37] = t
        mpra_type[k38] = t
        for tx, ok in (("RefSeq", gl.is_refseq(gene)), ("ENSEMBL", gl.is_ensembl(gene))):
            if not ok:
                continue
            for b, m, k in (("37", m37, k37), ("38", m38, k38)):
                coords[(tx, b)].setdefault(m[0], []).append(k)
                genes[(tx, b)].setdefault(k, []).append(gene)

    for chrom in CHROMS:
        for b in ("37", "38"):
            for tx in ("RefSeq", "ENSEMBL"):
                k = (tx, b)
                if chrom not in coords[k]:
                    continue
                lines = [f"{coord}\t{mpra_type[coord]}\t{';'.join(sorted(set(genes[k][coord])))}"
                         for coord in _sort_lines(sorted(set(coords[k][chrom])))]
                write_text_in_file("\n".join(lines), f"{out[k]}.tmp")

    for b in ("37", "38"):
        for tx in ("RefSeq", "ENSEMBL"):
            _sort_bed(f"{out[(tx, b)]}.tmp", out[(tx, b)], "checkMPRAfiles")

    if all(os.path.exists(p) for p in out.values()):
        _rm(prom_f)
        _rm(enh_f)
    if os.path.exists(necessary):
        c["MPRAann"] = 1


def check_regulatory_elements_files(config) -> None:
    """The bin/AnnotSV calls, in order."""
    check_promoter_file(config)
    check_ea_files(config)
    check_gh_files(config)
    check_mirtargetlink_files(config)
    check_abc_files(config)
    check_mpra_files(config)


# --------------------------------------------------------------------------
# regulatoryElementsAnnotation
# --------------------------------------------------------------------------
def sort_regulatory_elements_file(re_file: str) -> bool:
    """sortRegulatoryElementsFile."""
    env = dict(os.environ, LC_ALL="C")
    r = subprocess.run(["sort", "-c", "-s", "-k1,1", "-k2,2n", re_file], env=env,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if r.returncode == 0 and not r.stderr:
        print(f'\t...WARNING: "bedtools intersect -sorted" failed with {re_file}, the whole file is loaded into memory')
        return False
    if os.access(re_file, os.W_OK) and os.access(os.path.dirname(re_file) or ".", os.W_OK):
        print(f"\t...sorting {re_file} (not sorted by chromosome and start)")
        tmp = f"{re_file}.{os.getpid()}.tmp"
        with open(tmp, "wb") as out:
            r = subprocess.run(["sort", "-k1,1", "-k2,2n", re_file], env=env, stdout=out, stderr=subprocess.PIPE)
        if r.returncode == 0 and not r.stderr:
            os.replace(tmp, re_file)
            return True
        _rm(tmp)
    print(f'\t...WARNING: {re_file} is not sorted ("LC_ALL=C sort -k1,1 -k2,2n") and could not be sorted, '
          f"the whole file is loaded into memory")
    return False


@dataclass
class RegulatoryResult:
    # g_re: "chrom\tstart\tend" -> RE_gene text
    re: Dict[str, str] = field(default_factory=dict)
    # g_HI / g_TS of the regulated genes (ClinGen)
    hi: Dict[str, str] = field(default_factory=dict)
    ts: Dict[str, str] = field(default_factory=dict)
    all_regulated_genes: List[str] = field(default_factory=list)
    # g_PhenoGeniusCli / g_Exomiser (phenogenius.PhenoGenius, exomiser.Exomiser)
    phenogenius: object = None
    exomiser: object = None


class _Recorder:
    """recordRegulatedGenes state (g_genesReg, g_reDB, g_allRegulatedGenes...)."""

    def __init__(self, gl: GeneLists, clingen_hi, clingen_ts, skipped, select_on_read: bool):
        self.gl = gl
        self.hi, self.ts = clingen_hi, clingen_ts
        self.skipped = skipped
        self.select_on_read = select_on_read
        self.reportable: Dict[str, bool] = {}
        self.genes_reg: Dict[str, List[str]] = {}
        self.re_db: Dict[tuple, List[str]] = {}
        self.all_regulated: Dict[str, int] = {}

    def record(self, text: str) -> None:
        for L in _gets_lines(text):
            ls = _tsplit(L)
            sv = "\t".join(ls[0:3])
            db = ls[-2] if len(ls) >= 2 else ""
            for g in (ls[-1] if ls else "").split(";"):
                self.all_regulated[g] = 1
                if g in self.skipped:
                    continue
                if self.select_on_read:
                    if g not in self.reportable:
                        self.reportable[g] = (self.hi.get(g) == "3" or self.ts.get(g) == "3"
                                              or self.gl.is_morbid(g) or self.gl.is_candidate(g))
                    if not self.reportable[g]:
                        continue
                key = (sv, g)
                if key not in self.re_db:
                    self.genes_reg.setdefault(sv, []).append(g)
                    self.re_db[key] = [db]
                elif db not in self.re_db[key]:
                    self.re_db[key].append(db)


def _decode(data: bytes) -> str:
    """A Tcl pipe channel: system encoding, "auto" end of line translation."""
    text = data.decode("utf-8", errors="surrogateescape")
    return text.replace("\r\n", "\n").replace("\r", "\n")


def _strip_annotated_tsv(config, suffix: str) -> str:
    path = f"{config['outputDir']}/{config['outputFile']}"
    return re.sub(r"(.annotated)?.tsv$", lambda _m: suffix, path, count=1)


def regulatory_elements_annotation(config, all_genes_overlapped: List[str],
                                   phenogenius=None, exomiser=None) -> RegulatoryResult:
    """regulatoryElementsAnnotation.

    ``phenogenius`` / ``exomiser``: the ``phenogenius.PhenoGenius`` and
    ``exomiser.Exomiser`` objects, run here on every overlapped or regulated
    gene (with ``-hpo``), then used by OrganizeAnnotation.
    """
    c = config
    gl = _gene_lists(c)
    if phenogenius is None:
        from .phenogenius import PhenoGenius
        phenogenius = PhenoGenius(c)
    if exomiser is None:
        from .exomiser import Exomiser
        exomiser = Exomiser(c)
    pheno_on = str(c.get("PhenoGeniusCli", "0")) == "1"
    exo_on = str(c.get("Exomiser", "0")) == "1"
    reselect1 = str(c["REselect1"]) == "1"
    reselect2 = str(c["REselect2"]) == "1"
    rereport = str(c["REreport"]) == "1"

    # ClinGen HI/TS
    clingen_hi: Dict[str, str] = {}
    clingen_ts: Dict[str, str] = {}
    clingen_files = glob.glob(f"{c['annotationsDir']}/Annotations_{c['organism']}/Gene-based/ClinGen/"
                              f"*_ClinGenAnnotations.tsv*")
    if clingen_files:
        for L in _eof_loop_lines(_read(clingen_files[-1])):
            ls = _tsplit(L)
            gene, hi, ts = _li(ls, 0), _li(ls, 1), _li(ls, 2)
            if hi != "Not yet evaluated":
                clingen_hi[gene] = hi
            if ts != "Not yet evaluated":
                clingen_ts[gene] = ts

    select_on_read = reselect1 and not pheno_on and not exo_on
    skipped = set(all_genes_overlapped) if reselect2 else set()
    rec = _Recorder(gl, clingen_hi, clingen_ts, skipped, select_on_read)

    reg_dir = f"{_reg_dir(c)}/{c['genomeBuild']}"
    tx, build = c["tx"], c["genomeBuild"]
    re_files = [
        f"{reg_dir}/EA_{tx}_{build}.sorted.bed",
        f"{reg_dir}/GH_{tx}_{build}.sorted.bed",
        f"{reg_dir}/promoter_{c['promoterSize']}bp_{tx}_{build}.sorted.bed",
        f"{reg_dir}/miRTargetLink_{tx}_{build}.sorted.bed",
        f"{reg_dir}/ABC_{tx}_{build}.sorted.bed",
        f"{reg_dir}/MPRA_{tx}_{build}.sorted.bed",
    ]
    intersect_file = _strip_annotated_tsv(c, ".SV_RE_intersect.tmp")
    coords_file = _strip_annotated_tsv(c, ".SV_RE_coordinates.tmp.bed")
    stderr_file = _strip_annotated_tsv(c, ".SV_RE_intersect.stderr.tmp")
    _rm(intersect_file)
    if rereport:
        sv_bed = c["bedFile"]
    else:
        seen = set()
        with open(c["bedFile"], "r", errors="surrogateescape") as fin, \
                open(coords_file, "w", errors="surrogateescape") as fout:
            for L in _gets_lines(fin.read()):
                coord = "\t".join(L.split("\t")[0:3])
                if coord in seen:
                    continue
                seen.add(coord)
                fout.write(coord + "\n")
        sv_bed = coords_file

    bedtools = c.get("bedtools", "bedtools")
    message = ""
    for re_file in re_files:
        if not os.path.exists(re_file):
            continue
        ok = False
        if rereport:
            size = os.path.getsize(intersect_file) if os.path.exists(intersect_file) else 0
            for opt in (["-sorted"], []):
                if not opt and sort_regulatory_elements_file(re_file):
                    opt = ["-sorted"]
                with open(intersect_file, "ab") as out:
                    r = subprocess.run([bedtools, "intersect", *opt, "-a", sv_bed, "-b", re_file, "-wa", "-wb"],
                                       stdout=out, stderr=subprocess.PIPE)
                if r.returncode == 0 and not r.stderr:
                    ok = True
                    break
                message = r.stderr.decode(errors="replace").rstrip("\n") or "child process exited abnormally"
                if os.path.exists(intersect_file):
                    with open(intersect_file, "r+b") as f:
                        f.truncate(size)
        else:
            for opt in (["-sorted"], []):
                if not opt and sort_regulatory_elements_file(re_file):
                    opt = ["-sorted"]
                _rm(stderr_file)
                with open(stderr_file, "wb") as err:
                    r = subprocess.run([bedtools, "intersect", *opt, "-a", sv_bed, "-b", re_file, "-wa", "-wb"],
                                       stdout=subprocess.PIPE, stderr=err)
                rec.record(_decode(r.stdout))
                if r.returncode == 0 and os.path.getsize(stderr_file) == 0:
                    ok = True
                    break
                message = "child process exited abnormally"
            if not ok:
                try:
                    with open(stderr_file, "r", errors="replace") as f:
                        message += "\n" + f.read()
                except OSError:
                    pass
        if not ok:
            _rm(coords_file)
            raise ConfigExit(f"-- regulatoryElementsAnnotation --\n{bedtools} intersect -a {sv_bed} -b {re_file} -wa -wb\n"
                             f"{message}\nExit with error", 2)
    if rereport and os.path.exists(intersect_file):
        with open(intersect_file, "rb") as f:
            rec.record(_decode(f.read()))
    _rm(coords_file)
    _rm(stderr_file)

    all_regulated = sorted(rec.all_regulated)
    genes_reg = {sv: sorted(set(gs)) for sv, gs in rec.genes_reg.items()}
    re_db = {k: sorted(set(v)) for k, v in rec.re_db.items()}

    if os.path.exists(intersect_file):
        if rereport and os.path.getsize(intersect_file) != 0:
            report = _strip_annotated_tsv(c, ".SV_RE_intersect.report")
            _rm(report)
            os.rename(intersect_file, report)
        else:
            _rm(intersect_file)

    print("...searching for SV overlaps with a gene or a regulatory elements")
    print(f"\t...{len(all_genes_overlapped)} genes overlapped with an SV")
    print(f"\t...{len(all_regulated)} genes regulated by a regulatory element which is overlapped with an SV\n")
    sys.stdout.flush()

    all_genes = sorted(set(list(all_genes_overlapped) + all_regulated))
    if pheno_on:
        phenogenius.run(all_genes)
    if exo_on and all_genes:
        exomiser.run(all_genes)
    # A tool that failed to run is turned off (g_AnnotSV(PhenoGeniusCli/Exomiser) set to 0)
    pheno_on = str(c.get("PhenoGeniusCli", "0")) == "1"
    exo_on = str(c.get("Exomiser", "0")) == "1"

    result = RegulatoryResult(all_regulated_genes=all_regulated, phenogenius=phenogenius, exomiser=exomiser)
    for g in all_regulated:
        if g in clingen_hi:
            result.hi[g] = clingen_hi[g]
        if g in clingen_ts:
            result.ts[g] = clingen_ts[g]

    overlapped = set(all_genes_overlapped)
    for sv, genes in genes_reg.items():
        priority: List[str] = []
        other: List[str] = []
        noann: List[str] = []
        for g in genes:
            if reselect2 and g in overlapped:
                continue
            hi = result.hi.get(g, "")
            ts = result.ts.get(g, "")
            pg = phenogenius.annotation(g, "specificity") if pheno_on else ""
            ex = exomiser.annotation(g, "score") if exo_on else ""
            ann: List[str] = []
            prio = False
            if reselect1:
                if hi == "3":
                    ann.append(f"HI={hi}")
                if ts == "3":
                    ann.append(f"TS={ts}")
                if pg == "A":
                    ann.append("PG=A")
                    prio = True
                if ex != "" and tcl.compare(ex, "0.7") > 0:
                    ann.append(f"EX={ex}")
                    prio = True
                if gl.is_morbid(g):
                    ann.append("morbid")
                if gl.is_candidate(g):
                    ann.append("candidate")
                if not ann:
                    continue
                ann.append("RE=" + "+".join(re_db[(sv, g)]))
            else:
                if hi != "":
                    ann.append(f"HI={hi}")
                if ts != "":
                    ann.append(f"TS={ts}")
                if pg != "":
                    ann.append("PG=A")
                if ex not in ("", "0.0000", "-1.0"):
                    ann.append(f"EX={ex}")
                if gl.is_morbid(g):
                    ann.append("morbid")
                if gl.is_candidate(g):
                    ann.append("candidate")
                ann.append("RE=" + "+".join(re_db[(sv, g)]))
            if ann:
                (priority if prio else other).append(f"{g} ({'/'.join(ann)})")
            else:
                noann.append(g)
        items = priority or other or noann
        if len(items) > 50:
            result.re[sv] = ";".join(items[:50]) + "..."
        else:
            result.re[sv] = ";".join(items)
    return result
