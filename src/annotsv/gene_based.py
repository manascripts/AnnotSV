"""Gene-based annotations: port of the OMIM, gene lists and gene table modules.

Ported Tcl sources:

* ``AnnotSV-ncbi.tcl``: checkNCBIandHGNC
* ``AnnotSV-omim.tcl``: memorizeGeneNameAlias, checkOMIMfile, checkMorbidfile,
  fromOMIMtoPhenotype, fromOMIMgeneToPhenotype, isMorbid
* ``AnnotSV-gencc.tcl``: checkGenCCgeneFile, updateGenCCgeneFile
* ``AnnotSV-haploinsufficiency.tcl``: checkHIfile
* ``AnnotSV-clingen.tcl``: checkClinGenFile
* ``AnnotSV-loeuf-pLI.tcl``: checkLOEUFfile
* ``AnnotSV-exac.tcl``: checkGeneIntoleranceFile, checkCNVintoleranceFile
* ``AnnotSV-candidategenes.tcl``: isCandidate
* ``AnnotSV-extann.tcl``: ExternalAnnotations
* the "Gene-based" parts of ``OrganizeAnnotation`` (``AnnotSV-write.tcl``):
  the header columns and the per-line text, for split and full lines.

The ``check*`` procs build the ``Gene-based/*/*.tsv.gz`` files from the
downloads, once; ``ExternalAnnotations`` then reads every gene-based TSV
file and the write step reports one value per output column.
"""

from __future__ import annotations

import glob
import gzip
import os
import re
import subprocess
import time
from decimal import Decimal
from typing import Dict, List, Optional, Tuple

from . import tcl
from .config import ConfigExit

# --------------------------------------------------------------------------
# Small Tcl helpers
# --------------------------------------------------------------------------


def _read(path: str) -> str:
    """``read -nonewline`` with Tcl's default "auto" end of line translation."""
    if path.endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8", errors="surrogateescape") as f:
            text = f.read()
    else:
        with open(path, "r", encoding="utf-8", errors="surrogateescape") as f:
            text = f.read()
    return text[:-1] if text.endswith("\n") else text


def _split_lines(text: str) -> List[str]:
    """``split $text "\\n"``: an empty text gives an empty list."""
    return text.split("\n") if text != "" else []


def _lines_from_file(path: str) -> List[str]:
    """LinesFromFile (``""`` gives no line; a missing file raises, as ``open`` does)."""
    if path == "":
        return []
    with open(path, "r", encoding="utf-8", errors="surrogateescape") as f:
        text = f.read()
    return _split_lines(text[:-1] if text.endswith("\n") else text)


def _lines_from_gz(paths) -> List[str]:
    """LinesFromGZFile on a glob result: ``gzip -cd`` of every file, one after the other."""
    if isinstance(paths, str):
        paths = [paths] if paths else []
    if not paths:
        return []
    text = ""
    for p in paths:
        with gzip.open(p, "rt", encoding="utf-8", errors="surrogateescape") as f:
            text += f.read()
    return _split_lines(text[:-1] if text.endswith("\n") else text)


def _tsplit(line: str) -> List[str]:
    """``split $L "\\t"``."""
    return line.split("\t") if line != "" else []


def _li(ls: List[str], i: int) -> str:
    return ls[i] if 0 <= i < len(ls) else ""


def _lsearch_exact(ls: List[str], value: str) -> int:
    return ls.index(value) if value in ls else -1


def _lsearch_regexp(ls: List[str], pattern: str) -> int:
    for i, item in enumerate(ls):
        if re.search(pattern, item):
            return i
    return -1


def _lsort_unique(items) -> List[str]:
    return sorted(set(items))


def _now() -> str:
    return time.strftime("%B %d %Y - %H:%M")


def _today() -> str:
    return time.strftime("%Y%m%d")


def _write_text(text: str, path: str) -> None:
    """WriteTextInFile (append)."""
    with open(path, "a", encoding="utf-8", errors="surrogateescape", newline="") as f:
        f.write(text + "\n")


def _replace_text(text: str, path: str) -> None:
    """ReplaceTextInFile."""
    with open(path, "w", encoding="utf-8", errors="surrogateescape", newline="") as f:
        f.write(text + "\n")


def _gzip(path: str, proc: str) -> None:
    """``catch {exec gzip $path} Message``."""
    r = subprocess.run(["gzip", path], capture_output=True, text=True)
    if r.returncode != 0 or r.stderr:
        print(f"-- {proc} --")
        print(f"gzip {path}")
        print(f"{r.stderr.rstrip(chr(10))}\n")


def _delete(path) -> None:
    """``file delete -force`` on a glob result."""
    for p in ([path] if isinstance(path, str) else path):
        if p and os.path.isfile(p):
            os.remove(p)


def _keep_only_last(files: List[str], what: str) -> None:
    print(f"Several {what} files exist:")
    print(tcl.format_list(files))
    print(f"Keep only one: {files[-1]}\n")
    for f in files[:-1]:
        os.replace(f, f + ".notused")


def _is_true(value) -> bool:
    return str(value).strip().lower() in ("1", "yes", "true", "on")


_DOUBLE_RE = re.compile(r"(-)?[0-9]+\.[0-9]+(e[-+][0-9]+)?")


def tcl_double_str(value: float) -> str:
    """Tcl 8.6 text of a double: the shortest round-trip digits, in fixed
    notation for decimal exponents -5 < e < 17, else ``d.ddde+X``."""
    if value != value:
        return "NaN"
    if value in (float("inf"), float("-inf")):
        return "Inf" if value > 0 else "-Inf"
    s = "-" if str(value).startswith("-") else ""
    if value == 0:
        return s + "0.0"
    d = Decimal(repr(value)).normalize()
    sign, dg, exp = d.as_tuple()
    digits = "".join(map(str, dg))
    e10 = len(digits) + exp - 1
    if -5 < e10 < 17:
        if exp >= 0:
            return s + digits + "0" * exp + ".0"
        point = len(digits) + exp
        if point > 0:
            return s + digits[:point] + "." + digits[point:]
        return s + "0." + "0" * (-point) + digits
    mant = digits[0] + ("." + digits[1:] if len(digits) > 1 else "")
    return f"{s}{mant}e{'+' if e10 >= 0 else '-'}{abs(e10)}"


def tcl_expr_value(value: str) -> str:
    """``expr {$value}`` on a plain value: numbers are normalised, other
    strings come back unchanged."""
    if re.fullmatch(r"[+-]?[0-9]+", value):
        return str(int(value))
    n = tcl.to_number(value) if value.strip() == value and value != "" else None
    if isinstance(n, float):
        return tcl_double_str(n)
    return value


def csv_split(line: str, sep: str = ",", delim: str = '"') -> List[str]:
    """tcllib ``::csv::split $line $sep`` (standard, non ``-alternate`` mode)."""
    if line.endswith(sep + delim + delim):
        line = line[: -3] + sep + "\0" + delim + delim + "\0"
    if line.startswith(delim + delim + sep):
        line = "\0" + delim + delim + "\0" + sep + line[3:]
    # Tcl: regsub -all {(^${delChar}|${delChar}$)} ... is in braces and never matches
    keys = [(sep + delim * 3, sep + "\0" + delim), (delim * 3 + sep, delim + "\0" + sep),
            (delim * 2, delim), (delim, "\0")]
    line = _string_map(line, keys)
    i = 0
    while True:
        a = line.find("\0", i)
        if a < 0:
            break
        b = line.find("\0", a + 1)
        if b < 0:
            break
        rng = line[a:b + 1]
        if sep in rng:
            line = line[:a] + rng.replace(sep, "\1") + line[b + 1:]
        i = b + 1
    line = _string_map(line, [(sep, "\0"), ("\1", sep), ("\0", "")])
    return line.split("\0") if line != "" else []


def _string_map(text: str, pairs: List[Tuple[str, str]]) -> str:
    """Tcl ``string map``: at each position, the first matching key wins."""
    res = []
    i, n = 0, len(text)
    while i < n:
        for k, v in pairs:
            if text.startswith(k, i):
                res.append(v)
                i += len(k)
                break
        else:
            res.append(text[i])
            i += 1
    return "".join(res)


# --------------------------------------------------------------------------
# Paths and per-run state
# --------------------------------------------------------------------------


def _ann_dir(c) -> str:
    return f"{c['annotationsDir']}/Annotations_{c['organism']}"


def _gene_based_dir(c) -> str:
    return f"{_ann_dir(c)}/Gene-based"


def _state(config) -> dict:
    st = getattr(config, "_gene_based_state", None)
    if st is None:
        st = {}
        try:
            setattr(config, "_gene_based_state", st)
        except AttributeError:
            pass
    return st


# --------------------------------------------------------------------------
# checkNCBIandHGNC
# --------------------------------------------------------------------------


def check_ncbi_and_hgnc(config) -> None:
    c = config
    d = f"{_gene_based_dir(c)}/NCBIandHGNCgeneID"
    if not re.search("Human", c["organism"]):
        c["hpo"] = ""
    elif not os.path.exists(f"{d}/results.txt"):
        print("\nWARNING: No Exomiser annotations available.")
        print(f"...{d}/results.txt doesn't exist")
        c["hpo"] = ""
    elif not os.path.exists(f"{d}/geneSymbol_NCBIandHGNCgeneID.tsv"):
        out = ["genes\tNCBI_gene_ID\tHGNC_gene_ID"]
        seen = set()
        for L in _lines_from_file(f"{d}/results.txt"):
            Ls = _tsplit(L)
            ncbi = _li(Ls, 4)
            if ncbi == "" or re.search("NCBI Gene ID", ncbi, re.I):
                continue
            hgnc = _li(Ls, 0)
            approved = _li(Ls, 1)
            if approved != "" and approved not in seen:
                seen.add(approved)
                out.append(f"{approved}\t{ncbi}\t{hgnc}")
            previous = _li(Ls, 2).replace(" ", "")
            for ps in previous.split(","):
                if ps != "" and ps not in seen:
                    seen.add(ps)
                    out.append(f"{ps}\t{ncbi}\t{hgnc}")
            # Tcl parity: the "alias" loop walks the previous symbols again
            for as_ in previous.split(","):
                if as_ != "" and as_ not in seen:
                    seen.add(as_)
                    out.append(f"{as_}\t{ncbi}\t{hgnc}")
        _write_text("\n".join(out), f"{d}/geneSymbol_NCBIandHGNCgeneID.tsv")


# --------------------------------------------------------------------------
# memorizeGeneNameAlias
# --------------------------------------------------------------------------


def _overlaps(start: str, end: str, true_start: str, true_end: str) -> bool:
    return ((tcl.compare(start, true_start) < 0 and tcl.compare(end, true_start) > 0)
            or (tcl.compare(start, true_start) > 0 and tcl.compare(start, true_end) < 0))


def memorize_gene_name_alias(config) -> None:
    """Fill ``g_L_coord`` and ``g_alias`` (kept in the config's gene-based state)."""
    c = config
    st = _state(c)
    if c["organism"] != "Human":
        return
    for build in ("CHM13", "GRCh38", "GRCh37"):
        coord_file = f"{_ann_dir(c)}/Genes/{build}/genes.{c['tx']}.sorted.bed"
        if os.path.exists(coord_file):
            break
    coords: Dict[str, List[str]] = {}
    for L in _lines_from_file(coord_file):
        Ls = tcl.split_list(L)
        coords.setdefault(_li(Ls, 4), []).append(f"{_li(Ls, 0)} {_li(Ls, 1)} {_li(Ls, 2)}")
    coords = {g: _lsort_unique(v) for g, v in coords.items()}

    ncbi_for: Dict[str, str] = {}
    genes_for: Dict[str, List[str]] = {}
    for L in _lines_from_file(f"{_gene_based_dir(c)}/NCBIandHGNCgeneID/geneSymbol_NCBIandHGNCgeneID.tsv"):
        Ls = tcl.split_list(L)
        gene, ncbi = _li(Ls, 0), _li(Ls, 1)
        ncbi_for[gene] = ncbi
        genes_for.setdefault(ncbi, []).append(gene)

    alias: Dict[str, List[str]] = st.get("alias", {})
    for gene in coords:
        if gene not in ncbi_for:
            continue
        for coord_true in coords[gene]:
            t = tcl.split_list(coord_true)
            t_chrom, t_start, t_end = _li(t, 0), _li(t, 1), _li(t, 2)
            for alias_gene in genes_for[ncbi_for[gene]]:
                for coord in coords.get(alias_gene, []):
                    a = tcl.split_list(coord)
                    if t_chrom != _li(a, 0):
                        continue
                    if _overlaps(_li(a, 1), _li(a, 2), t_start, t_end):
                        alias.setdefault(gene, []).append(alias_gene)
                        break
    st["alias"] = {g: _lsort_unique(v) for g, v in alias.items()}
    st["coord"] = coords


# --------------------------------------------------------------------------
# checkOMIMfile
# --------------------------------------------------------------------------

_INHERITANCE = [("Autosomal dominant", "AD"), ("Autosomal recessive", "AR"), ("X-linked dominant", "XLD"),
                ("X-linked recessive", "XLR"), ("Y-linked dominant", "YLD"), ("Y-linked recessive", "YLR"),
                ("X-linked", "XL"), ("Y-linked", "YL")]
_PHENO_RE = re.compile(r" *(.*\([0-9]+\)),? *(.*)", re.S)


def check_omim_file(config) -> None:
    c = config
    st = _state(c)
    alias, coords = st.get("alias", {}), st.get("coord", {})
    omim_dir = f"{_gene_based_dir(c)}/OMIM"
    downloaded = glob.glob(f"{omim_dir}/genemap2.txt")
    f1 = glob.glob(f"{omim_dir}/*_OMIM-1-annotations.tsv.gz")
    f2 = glob.glob(f"{omim_dir}/*_OMIM-2-annotations.tsv.gz")
    if not downloaded:
        return
    if len(f1) > 1:
        _keep_only_last(f1, "OMIM")
        f1 = f1[-1:]
    if len(f2) > 1:
        _keep_only_last(f2, "OMIM")
        f2 = f2[-1:]
    if f1 and f2:
        return

    out1 = f"{omim_dir}/{_today()}_OMIM-1-annotations.tsv"
    out2 = f"{omim_dir}/{_today()}_OMIM-2-annotations.tsv"
    print(f"\t...creation of {os.path.basename(out1)}.gz and {os.path.basename(out2)} ({_now()})")
    print(f"\t   (in {omim_dir})")
    _replace_text("genes\tOMIM_ID", out1)
    _replace_text("genes\tOMIM_phenotype\tOMIM_inheritance", out2)

    idx: Dict[str, int] = {}
    pheno: Dict[str, List[str]] = {}
    inher: Dict[str, List[str]] = {}
    mims: Dict[str, List[str]] = {}
    l_genes: List[str] = []
    names = [("i_chrom", "# Chromosome", "# Chromosome"),
             ("i_start", "Genomic Position Start", "Genomic Position Start"),
             ("i_end", "Genomic Position End", "Genomic Position End"),
             ("i_gene1", "Gene/Locus And Other Related Symbols", "Gene/Locus And Other Related Symbols"),
             ("i_gene2", "Approved Gene Symbol", "Approved Symbol"),
             ("i_mim", "MIM Number", "Mim Number"),
             ("i_pheno", "Phenotypes", "Phenotypes")]
    for L in _lines_from_file(downloaded[0]):
        Ls = L.split("\t")
        if re.match("# Chromosome\tGenomic", L):
            for key, col, label in names:
                idx[key] = _lsearch_exact(Ls, col)
                if idx[key] == -1:
                    sep = "\n" if key == "i_pheno" else ":\n"
                    raise ConfigExit(f"{downloaded[0]}{sep}Bad header line syntax. \"{label}\" column not "
                                     "found - Exit with error", 2)
            continue
        if L.startswith("#"):
            continue
        if not idx:
            raise ConfigExit("can't read \"i_chrom\": no such variable", 1)
        true_chrom = _li(Ls, idx["i_chrom"]).replace("chr", "", 1)
        true_start, true_end = _li(Ls, idx["i_start"]), _li(Ls, idx["i_end"])
        gene1 = _li(Ls, idx["i_gene1"]).replace(" ", "")
        gene2 = _li(Ls, idx["i_gene2"])
        tmp = gene1.split(",") if gene1 != "" else []
        if gene2 != "":
            tmp.append(gene2)
        all_genes: List[str] = []
        for g in _lsort_unique(tmp):
            if g == "":
                continue
            for coord in coords.get(g, []):
                a = tcl.split_list(coord)
                if true_chrom != _li(a, 0):
                    continue
                if _overlaps(_li(a, 1), _li(a, 2), true_start, true_end):
                    all_genes.append(g)
                    all_genes.extend(alias.get(g, []))
                    break
        if not all_genes:
            continue
        all_genes = _lsort_unique(all_genes)

        mim = _li(Ls, idx["i_mim"])
        p_tmp = _li(Ls, idx["i_pheno"])
        for long, short in _INHERITANCE:
            p_tmp = re.sub(re.escape(long), short, p_tmp, flags=re.I)
        for p in p_tmp.split(";"):
            m = _PHENO_RE.search(p)
            if not m:
                continue
            p1, i1 = m.group(1), m.group(2)
            if i1 != "":
                i1 = re.sub(" *, +", ",", i1)
                for g in all_genes:
                    g = g.replace(" ", "")
                    if g == "":
                        continue
                    pheno.setdefault(g, []).append(f"{p1} {i1}")
                    inher.setdefault(g, []).extend(i1.split(";"))
            else:
                for g in all_genes:
                    g = g.replace(" ", "")
                    if g == "":
                        continue
                    pheno.setdefault(g, []).append(p1)
            for g in all_genes:
                g = g.replace(" ", "")
                if g == "":
                    continue
                l_genes.append(g)
                mims.setdefault(g, []).append(mim)

    for g in _lsort_unique(l_genes):
        ph = ";".join(_lsort_unique(pheno[g]))
        inh = ";".join(_lsort_unique(inher.get(g, [])))
        if re.fullmatch(";+", inh):
            inh = ""
        _write_text(f"{g}\t{';'.join(_lsort_unique(mims[g]))}", out1)
        _write_text(f"{g}\t{ph}\t{inh}", out2)
    _gzip(out1, "checkOMIMfile")
    _gzip(out2, "checkOMIMfile")
    _delete(downloaded)


# --------------------------------------------------------------------------
# checkMorbidfile
# --------------------------------------------------------------------------


def check_morbid_file(config) -> None:
    c = config
    st = _state(c)
    alias = st.get("alias", {})
    omim_dir = f"{_gene_based_dir(c)}/OMIM"
    downloaded = glob.glob(f"{omim_dir}/morbidmap.txt")
    formatted = glob.glob(f"{omim_dir}/*_morbid.tsv.gz")
    cand_formatted = glob.glob(f"{omim_dir}/*_morbidCandidate.tsv.gz")
    if not downloaded and not formatted:
        return
    if len(formatted) > 1:
        _keep_only_last(formatted, "Morbid Genes")
        return
    if len(cand_formatted) > 1:
        _keep_only_last(cand_formatted, "Morbid Genes candidate")
        return

    if not formatted:
        out = f"{omim_dir}/{_today()}_morbid.tsv"
        out_c = f"{omim_dir}/{_today()}_morbidCandidate.tsv"
        print(f"\t...creation of {os.path.basename(out)}.gz and {os.path.basename(out_c)} ({_now()})")
        print(f"\t   (in {omim_dir})")
        _replace_text("genes\tOMIM_morbid", out)
        _replace_text("genes\tOMIM_morbid_candidate", out_c)
        morbid: List[str] = []
        cand: List[str] = []
        i_gene = i_pheno = None
        for L in _lines_from_file(downloaded[0] if downloaded else ""):
            Ls = L.split("\t")
            if re.match("# Phenotype\t", L):
                i_gene = _lsearch_exact(Ls, "Gene/Locus And Other Related Symbols")
                if i_gene == -1:
                    raise ConfigExit("Bad header line syntax. \"Gene/Locus And Other Related Symbols\" column "
                                     "not found - Exit with error", 2)
                i_pheno = _lsearch_exact(Ls, "# Phenotype")
                if i_pheno == -1:
                    raise ConfigExit("Bad header line syntax. Phenotype column not found - Exit with error", 2)
                continue
            if L.startswith("#"):
                continue
            if i_gene is None:
                raise ConfigExit("can't read \"i_pheno\": no such variable", 1)
            pheno = _li(Ls, i_pheno)
            if re.match(r"\[.+\]", pheno, re.S):
                continue
            if re.search(r"\(1\)|\(2\)", pheno):
                continue
            gene = _li(Ls, i_gene)
            for g in _lsort_unique(gene.split(",") if gene != "" else []):
                g = g.replace(" ", "")
                if g == "":
                    continue
                target = cand if re.match(r"\{.+\}|\?", pheno, re.S) else morbid
                target.extend(alias[g] if g in alias else [g])
        for g in _lsort_unique(morbid):
            _write_text(f"{g}\tyes", out)
        _gzip(out, "checkMorbidfile")
        for g in _lsort_unique(cand):
            _write_text(f"{g}\tyes", out_c)
        _gzip(out_c, "checkMorbidfile")
        _delete(downloaded)

    # Tcl: "unset g_alias; unset g_L_coord"
    st.pop("alias", None)
    st.pop("coord", None)


# --------------------------------------------------------------------------
# fromOMIMtoPhenotype / fromOMIMgeneToPhenotype / isMorbid / isCandidate
# --------------------------------------------------------------------------


def _load_omim_phenotypes(config) -> None:
    st = _state(config)
    if "phen" in st:
        return
    omim_dir = f"{_gene_based_dir(config)}/OMIM"
    genes_of: Dict[str, List[str]] = {}
    for L in _lines_from_gz(glob.glob(f"{omim_dir}/*_OMIM-1-annotations.tsv.gz")):
        Ls = _tsplit(L)
        gene = _li(Ls, 0)
        for omim in (_li(Ls, 1).split(";") if _li(Ls, 1) != "" else []):
            genes_of.setdefault(omim, []).append(gene)
    phen: Dict[str, str] = {}
    for L in _lines_from_gz(glob.glob(f"{omim_dir}/*_OMIM-2-annotations.tsv.gz")):
        Ls = _tsplit(L)
        p = re.sub(r"\{|\}|\[|\]", "", _li(Ls, 1))
        phen[_li(Ls, 0)] = p.replace("/", ";")
    st["phen"], st["omim_genes"] = phen, genes_of


def from_omim_to_phenotype(config, omim: str) -> str:
    _load_omim_phenotypes(config)
    st = _state(config)
    found = [st["phen"][g] for g in st["omim_genes"].get(omim, []) if g in st["phen"]]
    return ";".join(found)


def from_omim_gene_to_phenotype(config, gene: str) -> str:
    _load_omim_phenotypes(config)
    return _state(config)["phen"].get(gene, "")


def is_morbid(config, gene: str) -> bool:
    """isMorbid."""
    st = _state(config)
    if "morbid" not in st:
        files = glob.glob(f"{_gene_based_dir(config)}/OMIM/*_morbid.tsv.gz")
        morbid = set()
        if files:
            text = ""
            for p in files:
                with gzip.open(p, "rt", encoding="utf-8", errors="surrogateescape") as f:
                    text += f.read()
            # Tcl parity: "while {![eof $f]}" also reads the empty line after the last
            # newline, so isMorbid "" is true
            for L in text.split("\n"):
                morbid.add(tcl.lindex(L, 0))
        st["morbid"] = morbid
    return gene in st["morbid"]


def is_candidate(config, gene: str) -> bool:
    """isCandidate."""
    st = _state(config)
    if "candidate" not in st:
        cand = set()
        path = config.get("candidateGenesFile", "") if hasattr(config, "get") else config["candidateGenesFile"]
        if path:
            with open(path, "r", encoding="utf-8", errors="surrogateescape") as f:
                text = f.read()
            for L in text.split("\n"):
                cand.update(tcl.split_chars(L, " |\n|\t"))
        st["candidate"] = cand
    return gene in st["candidate"]


# --------------------------------------------------------------------------
# checkHIfile, checkClinGenFile, checkLOEUFfile
# --------------------------------------------------------------------------


def check_hi_file(config) -> None:
    c = config
    d = f"{_gene_based_dir(c)}/DDD"
    downloaded = glob.glob(f"{d}/HI_Predictions_Version3.bed.gz")
    formatted = glob.glob(f"{d}/*_HI.tsv.gz")
    if not downloaded and not formatted:
        return
    if len(formatted) > 1:
        _keep_only_last(formatted, "Haploinsufficiency")
        return
    if formatted:
        return
    out = f"{d}/{_today()}_HI.tsv"
    print(f"\t...HI configuration ({_now()})")
    print(f"\t\t...creation of {out}.gz")
    print("\t\t   (done only once during the first HI annotation)")
    text = ["genes\tDDD_HI_percent"]
    for L in _lines_from_gz(downloaded):
        if L.startswith("track name"):
            continue
        infos = _li(_tsplit(L), 3)
        infos = infos.split("|") if infos != "" else []
        percent = _li(infos, len(infos) - 1).replace("%", "", 1)
        text.append(f"{_li(infos, 0)}\t{percent}")
    _write_text("\n".join(text), out)
    _gzip(out, "checkHIfile")
    _delete(downloaded)


def check_clingen_file(config) -> None:
    c = config
    d = f"{_gene_based_dir(c)}/ClinGen"
    downloaded = glob.glob(f"{d}/ClinGen_gene_curation_list_*.tsv")
    formatted = glob.glob(f"{d}/*_ClinGenAnnotations.tsv.gz")
    if not downloaded and not formatted:
        return
    if len(formatted) > 1:
        _keep_only_last(formatted, "ClinGen")
        return
    if formatted:
        return
    out = f"{d}/{_today()}_ClinGenAnnotations.tsv"
    print(f"...creation of {out}.gz ({_now()})\n")
    _replace_text("genes\tHI\tTS", out)
    idx = None
    seen = set()
    text: List[str] = []
    for L in _lines_from_file(downloaded[0] if downloaded else ""):
        Ls = L.split("\t")
        if L.startswith("#Gene Symbol"):
            idx = []
            for col, label in (("#Gene Symbol", "Gene Symbol"), ("Haploinsufficiency Score", None),
                               ("Triplosensitivity Score", None)):
                i = _lsearch_exact(Ls, col)
                if i == -1:
                    raise ConfigExit(f"Bad header line syntax. {label or col} column not found - Exit with error", 2)
                idx.append(i)
            continue
        if L.startswith("#"):
            continue
        if idx is None:
            raise ConfigExit("can't read \"i_gene\": no such variable", 1)
        gene, hi, ts = (_li(Ls, i) for i in idx)
        if gene in seen:
            continue
        seen.add(gene)
        hi = "" if hi == "Not yet evaluated" else hi
        ts = "" if ts == "Not yet evaluated" else ts
        text.append(f"{gene}\t{hi}\t{ts}")
    _write_text("\n".join(text), out)
    _gzip(out, "checkClinGenFile")
    _delete(downloaded)


def check_loeuf_file(config) -> None:
    c = config
    d = f"{_gene_based_dir(c)}/gnomAD"
    downloaded = glob.glob(f"{d}/gnomad.*.lof_metrics.by_gene.txt")
    formatted = glob.glob(f"{d}/*_gnomAD.LOEUF.pLI.annotations.tsv.gz")
    if not downloaded and not formatted:
        return
    if len(formatted) > 1:
        _keep_only_last(formatted, "LOEUF")
        return
    if formatted:
        return
    out = f"{d}/{_today()}_gnomAD.LOEUF.pLI.annotations.tsv"
    print(f"\t...creation of {out}.gz ({_now()})")
    print("\t   (done only once during the first LOEUF annotation)")
    _replace_text("genes\tLOEUF_bin\tpLI_gnomAD\tpLI_ExAC", out)
    loeuf: Dict[str, List[str]] = {}
    pli: Dict[str, List[str]] = {}
    pli_exac: Dict[str, List[str]] = {}
    idx = None
    path = downloaded[0]
    for L in _lines_from_file(path):
        Ls = L.split("\t")
        if L.startswith("gene"):
            idx = [0]
            for pat, name in (("^oe_lof_upper_bin", "oe_lof_upper_bin"), ("^pLI", "pLI"), ("^exac_pLI", "exac_pLI")):
                i = _lsearch_regexp(Ls, pat)
                if i == -1:
                    raise ConfigExit(f"Bad syntax into {path}.\n{name} field not found - Exit with error", 2)
                idx.append(i)
            continue
        if idx is None:
            raise ConfigExit("can't read \"i_gene\": no such variable", 1)
        gene, lo, pl, pe = (_li(Ls, i) for i in idx)
        if lo != "NA":
            loeuf.setdefault(gene, []).append(lo)
        if pl != "NA":
            pli.setdefault(gene, []).append(pl)
        if pe != "NA":
            pli_exac.setdefault(gene, []).append(pe)
    text = []
    for g in _lsort_unique(list(loeuf) + list(pli) + list(pli_exac)):
        # Tcl parity: plain (string) lsort, so "10" < "9"
        vals = [sorted(t[g])[0] if g in t else "" for t in (loeuf, pli, pli_exac)]
        text.append(f"{g}\t" + "\t".join(vals))
    _write_text("\n".join(text), out)
    _gzip(out, "checkLOEUFfile")
    _delete(downloaded)


# --------------------------------------------------------------------------
# checkGenCCgeneFile / updateGenCCgeneFile
# --------------------------------------------------------------------------


def check_gencc_gene_file(config, gene_lists) -> None:
    c = config
    d = f"{_gene_based_dir(c)}/GenCC"
    downloaded = glob.glob(f"{d}/submissions-export-tsv")
    formatted = glob.glob(f"{d}/*_GenCC.tsv") + glob.glob(f"{d}/*_GenCC.tsv.gz")
    if not downloaded and not formatted:
        return
    if len(formatted) > 1:
        _keep_only_last(formatted, "GenCC")
        return
    if not formatted:
        _update_gencc_gene_file(c, downloaded, gene_lists)


def _gencc_moi(moi: str) -> str:
    def has(p):
        return re.search(p, moi, re.I) is not None
    if has("Autosomal dominant"):
        if has("with maternal imprinting"):
            return "ADm"
        if has("with paternal imprinting"):
            return "ADp"
        return "AD"
    if has("Autosomal recessive"):
        return "AR"
    if has("Digenic"):
        return "2G"
    if has("Mitochondrial"):
        return "MT"
    if has("Semidominant"):
        return "sD"
    if has("Somatic mosaicism"):
        return "SOM"
    for x, p in (("X-linked", "XL"), ("Y-linked", "YL")):
        if has(x):
            if has("dominant"):
                return p + "D"
            if has("recessive"):
                return p + "R"
            return p
    return ""


def _update_gencc_gene_file(c, downloaded: List[str], gene_lists) -> None:
    d = f"{_gene_based_dir(c)}/GenCC"
    path = downloaded[0]
    out = f"{d}/{_today()}_GenCC.tsv"
    print("\t...GenCC genes configuration")
    print(f"\t\t...creation of {out}.gz")
    print("\t\t   (done only once during the first GenCC annotation)")
    text = ["genes\tGenCC_disease\tGenCC_moi\tGenCC_classification\tGenCC_pmid"]
    idx = None
    l_genes: List[str] = []
    disease: Dict[str, List[str]] = {}
    mois: Dict[str, List[str]] = {}
    classif: Dict[str, List[str]] = {}
    pmids: Dict[str, List[str]] = {}
    ntext = "field not found - Exit with error"
    for L in _lines_from_file(path):
        if L.startswith('"uuid'):
            Ls = csv_split(L, "\t")
            idx = []
            for name in ("gene_symbol", "disease_title", "moi_title", "classification_title", "submitted_as_pmids"):
                i = _lsearch_regexp(Ls, name)
                if i == -1:
                    raise ConfigExit(f"Bad syntax into {path}.\n{name} {ntext}", 2)
                idx.append(i)
            continue
        Ls = csv_split(L, "\t")
        if idx is None:
            raise ConfigExit("can't read \"i_gene\": no such variable", 1)
        gene, dis, moi, cl, pmid = (_li(Ls, i) for i in idx)
        l_genes.append(gene)
        if dis != "" and dis != '"':
            disease.setdefault(gene, []).append(dis)
        moi = _gencc_moi(moi)
        if moi != "":
            mois.setdefault(gene, []).append(moi)
            if re.search("incomplete penetrance|Reduced penetrance", L, re.I):
                mois[gene].append("IPVE")
        if cl != "" and cl != '"':
            classif.setdefault(gene, []).append(cl)
        pmid = pmid.replace('"', "")
        if pmid != "":
            pmid = pmid.replace(" ,", ";").replace(",", ";").replace(" ", ";")
            for p in pmid.split(";"):
                if p == "0" or p == "":
                    continue
                pmids.setdefault(gene, []).append(p)
    print("\t")
    for gene in _lsort_unique(l_genes):
        if not gene_lists.is_refseq(gene) and not gene_lists.is_ensembl(gene):
            continue
        vals = [";".join(_lsort_unique(t[gene])) if gene in t else "" for t in (disease, mois, classif, pmids)]
        text.append(f"{gene}\t" + "\t".join(vals))
    _write_text("\n".join(text), out)
    _gzip(out, "updateGenCCgeneFile")
    _delete(downloaded)


# --------------------------------------------------------------------------
# checkGeneIntoleranceFile / checkCNVintoleranceFile (ExAC)
# --------------------------------------------------------------------------


def check_gene_intolerance_file(config) -> None:
    c = config
    d = f"{_gene_based_dir(c)}/ExAC"
    downloaded = glob.glob(f"{d}/fordist_cleaned_nonpsych_z_pli_rec_null_data.txt")
    formatted = glob.glob(f"{d}/*_GeneIntolerance-Zscore.annotations.tsv.gz")
    if not downloaded and not formatted:
        return
    if len(formatted) > 1:
        _keep_only_last(formatted, "Gene Intolerance")
        return
    if formatted:
        return
    out = f"{d}/{_today()}_GeneIntolerance-Zscore.annotations.tsv"
    print(f"\t...GeneIntolerance configuration for Z scores annotation from ExAC ({_now()})")
    print(f"\t\t...creation of {out}.gz")
    print("\t\t   (done only once during the first GeneIntolerance annotation)")
    path = downloaded[0]
    text = ["genes\tsynZ_ExAC\tmisZ_ExAC"]
    idx = None
    for L in _lines_from_file(path):
        if L.startswith("transcript"):
            Ll = tcl.split_list(L)  # Tcl parity: lsearch on the line string
            idx = []
            for pat, name in (("^gene", "gene"), ("^syn_z", "syn_z"), ("^mis_z", "mis_z")):
                i = _lsearch_regexp(Ll, pat)
                if i == -1:
                    raise ConfigExit(f"Bad syntax into {path}.\n{name} field not found - Exit with error", 2)
                idx.append(i)
            continue
        if idx is None:
            raise ConfigExit("can't read \"i_gene\": no such variable", 1)
        Ls = _tsplit(L)
        text.append("\t".join(_li(Ls, i) for i in idx))
    _write_text("\n".join(text), out)
    _gzip(out, "checkGeneIntoleranceFile")
    _delete(downloaded)


def check_cnv_intolerance_file(config) -> None:
    c = config
    d = f"{_gene_based_dir(c)}/ExAC"
    downloaded = glob.glob(f"{d}/exac-final-cnv.gene.scores071316")
    formatted = glob.glob(f"{d}/*_ExAC.CNV-Zscore.annotations.tsv.gz")
    if not downloaded and not formatted:
        return
    if len(formatted) > 1:
        _keep_only_last(formatted, "CNV-Intolerant-Genes")
        return
    if formatted:
        return
    out = f"{d}/{_today()}_ExAC.CNV-Zscore.annotations.tsv"
    print(f"\t...ExAC CNV Intolerant Genes configuration ({_now()})")
    print(f"\t\t...creation of {out}.gz")
    print("\t\t   (done only once during the first CNV Intolerant Genes  annotation)")
    path = downloaded[0]
    text = ["genes\tdelZ_ExAC\tdupZ_ExAC\tcnvZ_ExAC"]
    idx = None
    scores: Dict[str, List[str]] = {}
    for L in _lines_from_file(path):
        Ls = L.split(" ") if L != "" else []
        if re.match("gene", tcl.format_list(Ls)):
            idx = []
            for pat, name in (("^gene_symbol", "gene_symbol"), ("^del.score", "del.score"),
                              ("^dup.score", "dup.score"), ("^cnv.score", "cnv.score")):
                i = _lsearch_regexp(Ls, pat)
                if i == -1:
                    raise ConfigExit(f"Bad syntax into {path}.\n{name} field not found - Exit with error", 2)
                idx.append(i)
            continue
        if idx is None:
            raise ConfigExit("can't read \"i_gene\": no such variable", 1)
        gene = _li(Ls, idx[0])
        vals = [_li(Ls, i) for i in idx[1:]]
        if gene not in scores:
            scores[gene] = vals
        else:
            scores[gene] = [v if tcl.compare(v, old) > 0 else old for v, old in zip(vals, scores[gene])]
    # Tcl writes these lines in "array names" (hash) order: a lookup-only cache, not compared
    for gene, vals in scores.items():
        text.append(f"{gene}\t" + "\t".join(vals))
    _write_text("\n".join(text), out)
    _gzip(out, "checkCNVintoleranceFile")
    _delete(downloaded)


# --------------------------------------------------------------------------
# All the gene-based check procs, in the bin/AnnotSV order
# --------------------------------------------------------------------------


def check_gene_based_files(config, gene_lists=None) -> None:
    """checkGenCCgeneFile ... checkLOEUFfile, as called from ``bin/AnnotSV``
    (the Exomiser and PhenoGenius checks are not ported yet)."""
    if gene_lists is None:
        from .regulatory_elements import _gene_lists
        gene_lists = _gene_lists(config)
    check_gencc_gene_file(config, gene_lists)
    check_ncbi_and_hgnc(config)
    memorize_gene_name_alias(config)
    check_omim_file(config)
    check_morbid_file(config)
    check_hi_file(config)
    check_clingen_file(config)
    check_loeuf_file(config)


def check_gene_intolerance_files(config) -> None:
    """checkGeneIntoleranceFile and checkCNVintoleranceFile (called after the benign SV checks)."""
    check_gene_intolerance_file(config)
    check_cnv_intolerance_file(config)


# --------------------------------------------------------------------------
# ExternalAnnotations and the gene-based output columns
# --------------------------------------------------------------------------


class ExternalAnnotations:
    """The ``g_ExtAnnotation`` array: every gene-based TSV file, loaded once."""

    def __init__(self, files: List[str]):
        self.files: List[str] = sorted(set(files))
        self.header: Dict[str, str] = {}
        self.data: Dict[str, Dict[str, str]] = {}
        self.display: List[str] = []
        for f in self.files:
            self._load(f)

    def _load(self, path: str) -> None:
        self.display.append(f"\t\t...{os.path.basename(path)}")
        lines = _lines_from_gz(path) if re.search(".gz$", path) else _lines_from_file(path)
        first = True
        nb_header = 0
        data: Dict[str, str] = {}
        self.header[path] = ""
        for L in lines:
            L = L.replace("{", "(").replace("}", ")")
            if not re.search("[a-z0-9]+", L, re.I) or L.startswith("#"):
                continue
            if first:
                first = False
                if not re.search("gene", L.split("\t")[0], re.I):
                    print(f"\t\tReading {path}:\n\tWARNING: Header is present ({L}) but does not contain a "
                          "\"gene\" named column.")
                else:
                    self.header[path] = L
                    nb_header = len(L.split("\t"))
                continue
            Ls = L.split("\t")
            gid = Ls[0].strip()
            if len(Ls) != nb_header:
                print(f"WARNING: {gid}, format (nb columns) is not good {len(tcl.split_list(L))} vs "
                      f"{nb_header} columns")
            if gid in data:
                print(f"WARNING: Reading Annotation file {os.path.basename(path)} and {gid} is seen multiple times.")
            data[gid] = "\t".join(Ls[1:])
        self.data[path] = data
        cols = _tsplit(self.header[path])[1:]
        self.display.append(f"\t\t({len(data)} gene identifiers and {len(_tsplit(self.header[path])) - 1} "
                            f"annotations columns: {', '.join(cols)})")

    def get(self, path: str, gene: str) -> str:
        return self.data.get(path, {}).get(gene, "")

    def nb_columns(self, path: str) -> int:
        return len(_tsplit(self.header[path])) - 1


class GeneBasedAnnotations:
    """The "Gene-based" header and text of OrganizeAnnotation."""

    def __init__(self, config, ext: Optional[ExternalAnnotations] = None):
        c = config
        self.config = c
        self.metrics_fr = c.get("metrics", "us") == "fr" if hasattr(c, "get") else c["metrics"] == "fr"
        self.enabled = _is_true(c["geneBasedAnn"]) if "geneBasedAnn" in c.values else bool(c["extann"])
        extann = c["extann"]
        if isinstance(extann, str):
            extann = tcl.split_list(extann)
        self.ext = ext if ext is not None else (ExternalAnnotations(extann) if extann else None)
        self.columns: List[str] = []
        self.indexes: List[int] = []
        if self.enabled and self.ext is not None:
            all_cols: List[str] = []
            for f in self.ext.files:
                h = re.sub("[^\t]+\t", "", self.ext.header[f], count=1)
                all_cols.extend(_tsplit(h))
            wanted = c["outputColHeader"]
            if isinstance(wanted, str):
                wanted = tcl.split_list(wanted)
            for j, col in enumerate(all_cols):
                if col in wanted:
                    self.columns.append(col)
                    self.indexes.append(j)
        if not self.indexes:
            self.enabled = False
            try:
                c["geneBasedAnn"] = "0"
            except Exception:  # pragma: no cover
                pass

    def header(self) -> str:
        return "".join(f"\t{col}" for col in self.columns) if self.enabled else ""

    def text(self, gene_name: str, mode: str) -> List[str]:
        """Values of the gene-based columns of one output line."""
        if not self.enabled:
            return []
        ext = self.ext
        values: List[str] = []
        if mode == "split":
            for f in ext.files:
                found = ext.get(f, gene_name)
                if found == "":
                    values.extend([""] * ext.nb_columns(f))
                elif self.metrics_fr:
                    for v in found.split("\t"):
                        if _DOUBLE_RE.fullmatch(v):
                            v = tcl_expr_value(v).replace(".", ",")
                        values.append(v)
                else:
                    values.extend(found.split("\t"))
        else:
            if gene_name == "":
                for f in ext.files:
                    values.extend([""] * ext.nb_columns(f))
            else:
                genes = gene_name.split(";")
                for f in ext.files:
                    values.extend(self._full_values(f, genes))
        return [values[i] if i < len(values) else "" for i in self.indexes]

    def _full_values(self, f: str, genes: List[str]) -> List[str]:
        from .general import merge_annotation
        ext = self.ext
        nb = ext.nb_columns(f)
        found = []
        for g in genes:
            a = ext.get(f, g)
            found.append(a if a != "" else "\t".join([""] * nb))
        merged = merge_annotation(found, nb)
        if merged == "":
            return [""]
        # Tcl parity: these regexps run on the whole path, not only on the file name
        is_clingen = re.search("ClinGenAnnotations.tsv", f) is not None
        is_omim1 = re.search("OMIM-1-annotations.tsv", f) is not None
        is_morbid_f = re.search("morbid", f) is not None
        is_loeuf = re.search("gnomAD.LOEUF.pLI.annotations.tsv", f) is not None
        out = []
        for column in merged.split("\t"):
            if column == "":
                out.append("")
                continue
            is_score = True
            mx = "-1000"
            for v in column.split(";"):
                if is_clingen:
                    # Tcl parity: "max" is reset for every gene, so the last gene wins unless one is 3
                    mx = ""
                    if v in ("3", "2", "1"):
                        if v == "3":
                            mx = "3"
                            break
                        if v == "2":
                            mx = "2"
                        elif mx != "2":
                            mx = "1"
                elif is_omim1:
                    if v != "":
                        mx = v if mx == "-1000" else f"{mx};{v}"
                elif is_morbid_f:
                    if "yes" in v:
                        mx = "yes"
                elif is_loeuf and re.fullmatch("[0-9]", v):
                    # Tcl parity: an integer pLI (0 or 1) also takes this "min" branch
                    if v != "" and (mx == "-1000" or tcl.compare(v, mx) < 0):
                        mx = v
                else:
                    if _DOUBLE_RE.fullmatch(v):
                        if tcl.compare(v, mx) > 0:
                            mx = v
                    else:
                        is_score = False
            if is_score:
                if self.metrics_fr:
                    mx = tcl_expr_value(mx).replace(".", ",")
                out.append(mx)
            else:
                out.append("")
        return out
