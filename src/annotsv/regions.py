"""Region annotations: port of the cytoband, gap, segdup, repeat, ENCODE
blacklist, TAD, user BED, GC content and COSMIC modules.

* ``AnnotSV-cytoband.tcl``        -> ``check_cytoband`` (annotated with ``UserBED``)
* ``AnnotSV-gap.tcl``, ``AnnotSV-segdup.tcl``, ``AnnotSV-repeat.tcl``,
  ``AnnotSV-encodeblacklist.tcl`` -> ``check_breakpoint_file`` and ``BreakpointAnnotation``
* ``AnnotSV-gccontent.tcl``       -> ``check_fasta_files`` and ``GCcontent``
* ``AnnotSV-tad.tcl``             -> ``check_tad_files`` and ``TADannotation``
* ``AnnotSV-userBED.tcl``         -> ``check_users_bed`` and ``UserBED``
* ``AnnotSV-cosmic.tcl``          -> ``check_cosmic_file``

``RegionAnnotations`` gathers the header and the per-line text that
``OrganizeAnnotation`` writes for them. The ``g_AnnotSV`` flags (``cytoband``,
``gapAnn``, ``segdupAnn``, ``repeatAnn``, ``ENCODEblacklistAnn``, ``tadAnn``,
``gcContentAnn``, ``COSMICann``) are kept in the config, as in Tcl.
"""

from __future__ import annotations

import glob
import os
import re
import shutil
import subprocess
import tarfile
import time
from typing import Dict, List, Optional

from . import tcl
from .config import ConfigExit
from .general import check_bed, first_line_from_file, remove_redundancy_without_sorting

KARYOTYPE = [str(i) for i in range(1, 23)] + ["X", "Y", "M", "MT"]


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def _ann_dir(c) -> str:
    return f"{c['annotationsDir']}/Annotations_{c['organism']}"


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


def _read_lines(path: str) -> List[str]:
    """The lines ``while {![eof $f]} {gets $f}`` reads (the last one may be "")."""
    with open(path, "r", encoding="utf-8", errors="surrogateescape", newline="") as f:
        return f.read().split("\n")


def _delete(path: str) -> None:
    if path and os.path.isdir(path) and not os.path.islink(path):
        shutil.rmtree(path)
    elif path and os.path.lexists(path):
        os.remove(path)


def _out_tmp(c, pattern: str, repl: str) -> str:
    """``regsub -nocase $pattern $g_AnnotSV(bedFile) $repl x`` then ``$outputDir/[file tail $x]``."""
    name = re.sub(pattern, lambda _m: repl, c["bedFile"], count=1, flags=re.I)
    return f"{c['outputDir']}/{os.path.basename(name)}"


def _unix_sort(c, src: str, dst: str, caller: str) -> None:
    """The ``sort.tmp.bash`` script: ``LC_ALL=C sort -k1,1 -k2,2n src > dst``."""
    env = dict(os.environ, LC_ALL="C")
    with open(dst, "wb") as out:
        r = subprocess.run(["sort", "-k1,1", "-k2,2n", src], stdout=out, stderr=subprocess.PIPE, env=env)
    if r.returncode or r.stderr:
        raise ConfigExit(f"-- {caller} --\nsort -k1,1 -k2,2n {src} > {dst}\n"
                         f"{r.stderr.decode(errors='replace').rstrip(chr(10))}\nExit with error", 2)


def _exec_ok(cmd: List[str], out_file: str) -> (bool, str):
    """``catch {exec ... > out_file}``: fails on a non-zero exit or any stderr output."""
    with open(out_file, "w") as out:
        r = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE)
    err = r.stderr.decode(errors="replace").rstrip("\n")
    return (r.returncode == 0 and not r.stderr), err


def _intersect(c, a: str, b: str, out_file: str, caller: str) -> None:
    """``bedtools intersect -sorted -a -b -wa -wb``, then without ``-sorted`` if it fails."""
    bt = c.get("bedtools", "bedtools")
    ok, _ = _exec_ok([bt, "intersect", "-sorted", "-a", a, "-b", b, "-wa", "-wb"], out_file)
    if ok:
        return
    ok, err = _exec_ok([bt, "intersect", "-a", a, "-b", b, "-wa", "-wb"], out_file)
    if not ok:
        raise ConfigExit(f"-- {caller} --\n{bt} intersect -sorted -a {a} -b {b} -wa -wb > {out_file}\n"
                         f"{err}\nExit with error", 2)


def _wanted(c) -> List[str]:
    cols = c["outputColHeader"]
    return tcl.split_list(cols) if isinstance(cols, str) else list(cols)


def _any_wanted(c, cols) -> bool:
    wanted = _wanted(c)
    return any(col in wanted for col in cols)


def _flag(c, key: str) -> bool:
    return str(c.get(key, "0")) == "1"


def _keep_only_last(files: List[str], what: str) -> None:
    print(f"Several {what} files exist:")
    print(tcl.format_list(files))
    print(f"Keep only one: {files[-1]}\n")
    for f in files[:-1]:
        os.replace(f, f + ".notused")


def number_of_annotation_col(config) -> Dict[str, int]:
    """``g_numberOfAnnotationCol``."""
    st = config.get("_numberOfAnnotationCol")
    if st is None:
        st = {}
        config["_numberOfAnnotationCol"] = st
    return st


# --------------------------------------------------------------------------
# Cytoband
# --------------------------------------------------------------------------
def cytoband_dir(c) -> str:
    return f"{_ann_dir(c)}/AnyOverlap/CytoBand/{c['genomeBuild']}"


def check_cytoband(config) -> None:
    """checkCytoband."""
    c = config
    d = cytoband_dir(c)
    build = c["genomeBuild"]
    bed = f"{d}/cytoBand_{build}.bed"
    formatted = f"{d}/cytoBand_{build}.formatted.bed"
    formatted_sorted = f"{d}/cytoBand_{build}.formatted.sorted.bed"
    header_file = f"{d}/cytoBand_{build}.header.tsv"
    if os.path.exists(bed):
        check_bed(bed, c["outputDir"], d)
        _write_text("#chrom\tstart\tend\tCytoBand", header_file)
        _unix_sort(c, formatted, formatted_sorted, "checkCytoband")
        _delete(formatted)
        _delete(bed)
    number_of_annotation_col(c)[formatted_sorted] = 1
    c["cytoband"] = "1" if os.path.exists(formatted_sorted) and os.path.exists(header_file) else "0"


# --------------------------------------------------------------------------
# Breakpoint annotations: Gap, SegDup, Repeat, ENCODE blacklist
# --------------------------------------------------------------------------
BREAKPOINT_SOURCES = {
    # name: (flag, output columns, "Several ... files exist" label, header of the formatted file)
    "Repeat": ("repeatAnn", ("Repeat_coord_left", "Repeat_type_left", "Repeat_coord_right", "Repeat_type_right"),
               "repeat", None),
    "SegDup": ("segdupAnn", ("SegDup_left", "SegDup_right"), "segmental duplication", "#Chrom\tStart\tEnd\tSegDup"),
    "Gap": ("gapAnn", ("Gap_left", "Gap_right"), "gap", "#Chrom\tStart\tEnd\tGap"),
    "ENCODEblacklist": ("ENCODEblacklistAnn",
                        ("ENCODE_blacklist_left", "ENCODE_blacklist_characteristics_left",
                         "ENCODE_blacklist_right", "ENCODE_blacklist_characteristics_right"),
                        "ENCODE blacklist", "#Chrom\tStart\tEnd\tENCODE_blacklist\tENCODE_blacklist_characteristics"),
}
_REPEAT_CHROMS = [str(i) for i in range(1, 25)] + ["X", "Y", "M", "MT"]


def _bp_dir(c, name: str) -> str:
    return f"{_ann_dir(c)}/BreakpointsAnnotations/{name}/{c['genomeBuild']}"


def _bp_formatted(c, name: str) -> List[str]:
    return glob.glob(f"{_bp_dir(c, name)}/*_{name}.sorted.bed")


def check_breakpoint_file(config, name: str) -> None:
    """checkRepeatFile, checkSegDupFile, checkGapFile and checkENCODEblacklistFile."""
    c = config
    flag, cols, label, header = BREAKPOINT_SOURCES[name]
    d = _bp_dir(c, name)
    downloaded = glob.glob(f"{d}/{name}.bed")
    formatted = _bp_formatted(c, name)
    if not downloaded and not formatted:
        c[flag] = "0"
        return
    c[flag] = "1"
    if not _any_wanted(c, cols):
        c[flag] = "0"
        return
    if len(formatted) > 1:
        _keep_only_last(formatted, label)
        return
    if formatted:
        return

    # The downloaded file exists but not the formatted one
    out = f"{d}/{_today()}_{name}.sorted.bed"
    if name == "Repeat":
        _delete(f"{out}.tmp")
        print(f"\t...Repeats configuration ({_now()})")
        print(f"\t\t...creation of {out}")
        print("\t\t   (done only once during the first Repeats annotation)")
    else:
        print(f"\t...creation of the \"{os.path.basename(out)}\" file ({_now()})\n\t   (done only once)")
        _replace_text(header, f"{out}.tmp")

    text: List[str] = []
    i = 0
    for line in _read_lines(downloaded[0]):
        if line == "":
            continue
        if name == "Repeat":
            line = re.sub("^chr", "", line, count=1)
            if tcl.lindex(line, 0) not in _REPEAT_CHROMS:
                continue
            text.append(line)
            # Tcl parity: checkRepeatFile never increments $i, so the lines are written all at once
        else:
            if line.startswith("#"):
                continue
            ls = line.split("\t")
            chrom = ls[0].replace("chr", "", 1)
            start = ls[1] if len(ls) > 1 else ""
            end = ls[2] if len(ls) > 2 else ""
            row = f"{chrom}\t{start}\t{end}\t{chrom}:{start}-{end}"
            if name == "ENCODEblacklist":
                row += f"\t{ls[3] if len(ls) > 3 else ''}"
            text.append(row)
            if i > 500000:
                _write_text("\n".join(text), f"{out}.tmp")
                text = []
                i = 0
            i += 1
    _write_text("\n".join(text), f"{out}.tmp")
    _unix_sort(c, f"{out}.tmp", out, f"check{name}File")
    _delete(f"{out}.tmp")
    _delete(downloaded[0])


def check_breakpoint_files(config) -> None:
    """The breakpoint checks of bin/AnnotSV (after checkFASTAfiles), in that order."""
    for name in ("Repeat", "SegDup", "Gap", "ENCODEblacklist"):
        check_breakpoint_file(config, name)


def breakpoints_file(config) -> str:
    return _out_tmp(config, r"(.formatted)?.bed$", ".breakpoints.bed")


def _write_breakpoints(c, path: str, size_of: Optional[Dict[str, int]]) -> None:
    """The breakpoints +/- 100 bp file. ``size_of`` (GC content only) clips
    the windows at the end of the chromosome."""
    lines: List[str] = []
    i = 0
    tmp = f"{path}.tmp"
    for line in _read_lines(c["bedFile"]):
        if line == "":
            continue
        ls = line.split("\t")
        chrom, start, end = ls[0], int(ls[1]), int(ls[2])
        s1, e1 = max(start - 100, 0), max(start + 100, 0)
        s2, e2 = max(end - 100, 0), max(end + 100, 0)
        if size_of is not None and chrom in size_of:
            size = size_of[chrom]
            if tcl.compare(str(e1), str(size)) > 0 and tcl.compare(str(s1), str(size)) < 0:
                e1 = size
            if tcl.compare(str(e2), str(size)) > 0 and tcl.compare(str(s2), str(size)) < 0:
                e2 = size
        lines.append(f"{chrom}\t{s1}\t{e1}\t{ls[1]}")
        lines.append(f"{chrom}\t{s2}\t{e2}\t{ls[2]}")
        i += 1
        if i > 100000:
            _write_text("\n".join(sorted(set(lines))), tmp)
            lines = []
            i = 0
    _write_text("\n".join(sorted(set(lines))), tmp)


class BreakpointAnnotation:
    """GapAnnotation, SegDupAnnotation, RepeatAnnotation and ENCODEblacklistAnnotation."""

    def __init__(self, config, name: str):
        self.config = config
        self.name = name
        self.coord: Optional[Dict[str, List[str]]] = None
        self.other: Dict[str, List[str]] = {}
        # Tcl parity: one value for Gap/SegDup, "coord\ttype" for Repeat and ENCODE blacklist
        self.two_cols = name in ("Repeat", "ENCODEblacklist")
        self.empty = "\t" if self.two_cols else ""

    def _load(self) -> None:
        c = self.config
        formatted = _bp_formatted(c, self.name)
        tmp_bp = breakpoints_file(c)
        if not os.path.exists(tmp_bp):
            _delete(f"{tmp_bp}.tmp")
            _write_breakpoints(c, tmp_bp, None)
            _unix_sort(c, f"{tmp_bp}.tmp", tmp_bp, f"{self.name}Annotation, sort")
            _delete(f"{tmp_bp}.tmp")
        tmp = _out_tmp(c, r"(.formatted)?.bed$", f".intersect.{self.name}.bed")
        _delete(tmp)
        # Tcl: "-b $FileFormatted" is the glob result (a single file, checked before)
        _intersect(c, tmp_bp, formatted[0] if formatted else "", tmp, f"{self.name}Annotation, intersect")
        self.coord = {}
        for line in _read_lines(tmp):
            if line == "":
                continue
            ls = line.split("\t")
            key = f"{ls[0]},{ls[3]}"

            def g(j):
                return ls[j] if j < len(ls) else ""
            if self.name == "Repeat":
                self.coord.setdefault(key, []).append(f"{ls[0]}:{g(5)}-{g(6)}")
                self.other.setdefault(key, []).append(g(7))
            elif self.name == "ENCODEblacklist":
                self.coord.setdefault(key, []).append(g(7))
                self.other.setdefault(key, []).append(g(8))
            else:
                self.coord.setdefault(key, []).append(g(7))
        _delete(tmp)

    def __call__(self, chrom: str, pos: str) -> str:
        if self.coord is None:
            self._load()
        key = f"{chrom},{pos}"
        if key not in self.coord:
            return self.empty
        text = ";".join(self.coord[key])
        if self.two_cols:
            text += "\t" + ";".join(self.other[key])
        return text


# --------------------------------------------------------------------------
# GC content
# --------------------------------------------------------------------------
def _gc_dir(c) -> str:
    return f"{_ann_dir(c)}/BreakpointsAnnotations/GCcontent"


def check_fasta_files(config) -> None:
    """checkFASTAfiles."""
    c = config
    build = c["genomeBuild"]
    downloaded = glob.glob(f"{_gc_dir(c)}/{build}/*chromFa.tar.gz")
    formatted = glob.glob(f"{_gc_dir(c)}/{build}/{build}_chromFa.fasta")
    if not downloaded and not formatted:
        c["gcContentAnn"] = "0"
        return
    c["gcContentAnn"] = "1"
    if not _any_wanted(c, ("GC_content_left", "GC_content_right")):
        c["gcContentAnn"] = "0"
        return
    if not formatted:
        update_fasta_files(c)


def _sed_chr(path: str) -> str:
    """``exec sed "s/chr//" $path``: the first "chr" of each line removed, one trailing newline dropped."""
    with open(path, "r", encoding="latin-1", newline="") as f:
        text = f.read()
    text = "\n".join(line.replace("chr", "", 1) for line in text.split("\n"))
    return text[:-1] if text.endswith("\n") else text


def update_fasta_files(config) -> None:
    """updateFASTAfiles."""
    c = config
    build = c["genomeBuild"]
    gc_dir = _gc_dir(c)
    downloaded = glob.glob(f"{gc_dir}/{build}/*chromFa.tar.gz")
    formatted = f"{gc_dir}/{build}/{build}_chromFa.fasta"
    sizes = f"{gc_dir}/{build}/{build}_chromFa.chromSizes"
    print(f"\t...GC content configuration ({_now()})")
    print(f"\t\t...creation of {formatted}")
    print("\t\t   (done only once during the first GC content annotation)")

    # Tcl: gzip -cd | ::tar::untar -chan
    for path in downloaded:
        with tarfile.open(path, "r:gz") as tar:
            tar.extractall(f"{gc_dir}/{build}", **({"filter": "tar"} if hasattr(tarfile, "data_filter") else {}))

    path_gb = None
    for chrom in KARYOTYPE:
        if build == "GRCh37" or build in ("mm9", "mm10", "mm39"):
            path_gb = build
        elif re.search("GRCh38|CHM13", build):
            path_gb = f"{build}/chroms"
        if path_gb is None:
            continue
        fa = f"{gc_dir}/{path_gb}/chr{chrom}.fa"
        if not os.path.exists(fa):
            continue
        with open(formatted, "a", encoding="latin-1", newline="") as f:
            f.write(_sed_chr(fa) + "\n")
        if not re.search("Y|M", chrom):
            with open(fa, "rb") as f:
                data = f.read()
            n_chars, n_lines = len(data), data.count(b"\n")
            _write_text(f"{chrom}\t{n_chars - n_lines - len('>chr' + chrom)}", sizes)

    if path_gb is not None:
        for tmp in glob.glob(f"{gc_dir}/{path_gb}/chr*.fa"):
            _delete(tmp)
    _delete(f"{gc_dir}/{build}/chroms")

    # index file generation
    test_bed = f"{gc_dir}/test.bed"
    _replace_text("1\t150000\t150200", test_bed)
    subprocess.run(["bedtools", "nuc", "-fi", formatted, "-bed", test_bed], capture_output=True)
    _delete(test_bed)
    for path in downloaded:
        _delete(path)


class GCcontent:
    """GCcontentAnnotation."""

    def __init__(self, config):
        self.config = config
        self.values: Optional[Dict[str, str]] = None

    def _load(self) -> None:
        c = self.config
        build = c["genomeBuild"]
        size_of: Dict[str, int] = {}
        sizes = f"{_gc_dir(c)}/{build}/{build}_chromFa.chromSizes"
        with open(sizes) as f:
            for line in f.read().rstrip("\n").split("\n"):
                items = tcl.split_list(line)
                if items:
                    size_of[items[0]] = int(items[1]) if len(items) > 1 else 0
        fasta = " ".join(glob.glob(f"{_gc_dir(c)}/{build}/{build}_chromFa.fasta"))

        tmp_bp = breakpoints_file(c)
        _delete(f"{tmp_bp}.tmp")
        _delete(tmp_bp)
        _write_breakpoints(c, tmp_bp, size_of)
        bt = c.get("bedtools", "bedtools")
        ok, err = _exec_ok([bt, "sort", "-i", f"{tmp_bp}.tmp"], tmp_bp)
        if not ok:
            print(f"-- GCcontentAnnotation, sort --\n{bt} sort -i {tmp_bp}.tmp > {tmp_bp}\n{err}")
        _delete(f"{tmp_bp}.tmp")

        nuc = _out_tmp(c, r"(.formatted)?.bed$", ".GCcontent.txt")
        _delete(nuc)
        ok, err = _exec_ok([bt, "nuc", "-fi", fasta, "-bed", tmp_bp], nuc)
        if not ok and err != "Warning: the index file is older than the FASTA file.":
            print(f"-- GCcontentAnnotation, nuc --\n{bt} nuc -fi {fasta} -bed {tmp_bp} > {nuc}\n{err}")

        self.values = {}
        i_gc = i_n = i_oth = -1
        fr = c.get("metrics", "us") == "fr"
        for line in _read_lines(nuc):
            if line == "":
                continue
            ls = line.split("\t")
            if "_pct_gc" in line:
                i_gc = next((j for j, v in enumerate(ls) if "_pct_gc" in v), -1)
                i_n = next((j for j, v in enumerate(ls) if "_num_N" in v), -1)
                i_oth = next((j for j, v in enumerate(ls) if "_num_oth" in v), -1)
                continue
            key = f"{ls[0]},{ls[3]}"
            gc = f"{float(ls[i_gc]):.3f}"
            num_n, num_oth = ls[i_n], ls[i_oth]
            if fr:
                gc = gc.replace(".", ",")
            if num_n != "0":
                gc += f" ({num_n} N, {num_oth} Oth)" if num_oth != "0" else f" ({num_n} N)"
            elif num_oth != "0":
                gc += f" ({num_oth} Oth)"
            self.values[key] = gc
        _delete(nuc)

    def __call__(self, chrom: str, pos: str) -> str:
        if self.values is None:
            self._load()
        return self.values.get(f"{chrom},{pos}", "")


# --------------------------------------------------------------------------
# TAD
# --------------------------------------------------------------------------
def _tad_dir(c) -> str:
    return f"{_ann_dir(c)}/FtIncludedInSV/TAD/{c['genomeBuild']}"


def _tad_formatted(c) -> List[str]:
    return glob.glob(f"{_tad_dir(c)}/*_boundariesTAD.sorted.bed")


def _tcl_int(value: str) -> int:
    return int(tcl.to_number(value))


def check_tad_files(config) -> None:
    """checkTADfiles."""
    c = config
    downloaded = glob.glob(f"{_tad_dir(c)}/ENC*.bed.gz")
    formatted = _tad_formatted(c)
    if not downloaded and not formatted:
        c["tadAnn"] = "0"
        return
    c["tadAnn"] = "1"
    if not _any_wanted(c, ("TAD_coordinate", "ENCODE_experiment")):
        c["tadAnn"] = "0"
        return
    if len(formatted) > 1:
        _keep_only_last(formatted, "TAD")
        return
    if formatted:
        return

    out = f"{_tad_dir(c)}/{_today()}_boundariesTAD.sorted.bed"
    print("\t...TAD configuration")
    print(f"\t\t...creation of {out}")
    print("\t\t   (done only once during the first boundaries TAD annotation)")
    experiment: Dict[str, List[str]] = {}
    tad: Dict[str, List[str]] = {}
    all_lines: Dict[str, List[str]] = {}
    for tad_file in glob.glob(f"{_tad_dir(c)}/ENC*.bed.gz"):
        # Tcl parity: regsub ".bed$" does not match "*.bed.gz", so the experiment keeps ".bed.gz"
        exp = re.sub(r".bed$", "", tad_file, count=1, flags=re.I)
        exp = re.sub(r".*/", "", exp, count=1)
        r = subprocess.run(["gzip", "-cd", tad_file], capture_output=True)
        content = r.stdout.decode("utf-8", errors="surrogateescape")
        if content.endswith("\n"):
            content = content[:-1]
        for line in content.split("\n"):
            if line.startswith("#"):
                continue
            if line == "":
                # Tcl parity: an empty line ends on "if {"" <= 0} continue"
                continue
            items = tcl.split_list(line)

            def li(j):
                return items[j] if j < len(items) else ""
            chrom = li(0).replace("chr", "", 1)
            tad_start, tad_end = li(1), li(2)
            if "___boundary" in line:
                b_start = tad_start
                b_end = str(_tcl_int(tad_start) + 40000)
            else:
                b_start = str(_tcl_int(tad_start) - 40000)
                b_end = tad_start
                tad_start = b_start
                tad_end = str(_tcl_int(tad_end) + 40000)
            if tcl.compare(b_start, "0") < 0:
                b_start = "0"
            if tcl.compare(b_end, b_start) <= 0:
                continue
            key = f"{chrom}\t{b_start}\t{b_end}"
            experiment.setdefault(key, []).append(exp)
            tad.setdefault(key, []).append(f"{chrom}:{tad_start}-{tad_end}")
            all_lines.setdefault(chrom, []).append(f"{b_start}\t{b_end}")

    for chrom in KARYOTYPE:
        if chrom not in all_lines:
            continue
        bounds_list = tcl.sort_on_element(tcl.sort_on_element(sorted(set(all_lines[chrom])), 1), 0)
        for bounds in bounds_list:
            key = f"{chrom}\t{bounds}"
            _write_text(f"{key}\t{','.join(sorted(set(experiment[key])))}\t{','.join(sorted(set(tad[key])))}",
                        f"{out}.tmp")
    _unix_sort(c, f"{out}.tmp", out, "checkTADfiles")
    _delete(f"{out}.tmp")
    for f in downloaded:
        _delete(f)


class TADannotation:
    """TADannotation: "TAD_coordinate\tENCODE_experiment", restricted to the requested columns."""

    def __init__(self, config):
        self.config = config
        self.indexes: List[int] = []
        self.columns: List[str] = []
        self.text: Optional[Dict[str, str]] = None
        self.empty = ""

    def header(self) -> str:
        """The "TAD header" part of OrganizeAnnotation (sets ``tadAnn_i``)."""
        c = self.config
        if not _flag(c, "tadAnn"):
            return ""
        wanted = _wanted(c)
        self.indexes, self.columns = [], []
        for j, col in enumerate(("TAD_coordinate", "ENCODE_experiment")):
            if col in wanted:
                self.columns.append(col)
                self.indexes.append(j)
        if not self.indexes:
            c["tadAnn"] = "0"
            return ""
        return "".join(f"\t{col}" for col in self.columns)

    def _load(self) -> None:
        c = self.config
        self.empty = "\t".join("" for _ in self.indexes)
        formatted = _tad_formatted(c)
        tmp = _out_tmp(c, r"(.formatted)?.bed$", ".intersect.tad")
        _delete(tmp)
        _intersect(c, c["bedFile"], formatted[0] if formatted else "", tmp, "TADannotation")
        experiments: Dict[str, List[str]] = {}
        coords: Dict[str, List[str]] = {}
        for line in _read_lines(tmp):
            if line == "":
                continue
            ls = line.split("\t")
            sv_start, sv_end = int(ls[1]), int(ls[2])
            b_start, b_end = int(ls[-4]), int(ls[-3])
            b_length = b_end - b_start
            sv_length = sv_end - sv_start
            if sv_length > 1:
                overlap_start = b_start if sv_start < b_start else sv_start
                overlap_end = sv_end if sv_end < b_end else b_end
                if (overlap_end - overlap_start) * 100.0 / b_length < 100:
                    continue
            key = f"{ls[0]},{ls[1]},{ls[2]}"
            experiments.setdefault(key, []).extend(ls[-2].split(","))
            coords.setdefault(key, []).extend(ls[-1].split(","))
        self.text = {}
        for key in experiments:
            full = [";".join(sorted(set(coords[key]))[0:21]), ":".join(sorted(set(experiments[key]))[0:21])]
            self.text[key] = "\t".join(full[i] for i in self.indexes)
        _delete(tmp)

    def __call__(self, chrom: str, start: str, end: str) -> str:
        if self.text is None:
            self._load()
        return self.text.get(f"{chrom},{start},{end}", self.empty)


# --------------------------------------------------------------------------
# COSMIC
# --------------------------------------------------------------------------
def check_cosmic_file(config) -> None:
    """checkCOSMICfile: build ``Users/<build>/FtIncludedInSV/CosmicCompleteCNA_<build>.bed``."""
    c = config
    build = c["genomeBuild"]
    downloaded = f"{_ann_dir(c)}/FtIncludedInSV/COSMIC/{build}/CosmicCompleteCNA.tsv.gz"
    users_dir = f"{_ann_dir(c)}/Users/{build}/FtIncludedInSV"
    formatted = f"{users_dir}/CosmicCompleteCNA_{build}.bed"
    if os.path.exists(formatted):
        c["COSMICann"] = "1"
        return
    if not os.path.exists(downloaded):
        return
    print(f"\t...COSMIC configuration ({_now()})")
    print(f"\t\t...reading {downloaded}")
    r = subprocess.run(["gzip", "-cd", downloaded], capture_output=True)
    lines = r.stdout.decode("utf-8", errors="surrogateescape").split("\n")
    to_write: Dict[str, List[str]] = {}
    header = None
    i = 0
    i_id = i_type = i_coord = i_chrom = i_start = i_end = -1
    chrom = start = end = None

    def lsearch_re(ls, pattern):
        return next((j for j, v in enumerate(ls) if re.search(pattern, v)), -1)

    def li(ls, j):
        return ls[j] if 0 <= j < len(ls) else ""
    for line in lines:
        ls = line.split("\t")
        i += 1
        if i % 10000000 == 0:
            print(f"\t{i} lines read")
        if re.search("^COSMIC_CNV_ID|CNV_ID", line):
            i_id = lsearch_re(ls, "^COSMIC_CNV_ID|CNV_ID")
            i_type = ls.index("MUT_TYPE") if "MUT_TYPE" in ls else -1
            i_coord = lsearch_re(ls, "Chromosome:G_Start")   # Old raw data version (before 2023)
            i_chrom = lsearch_re(ls, "CHROMOSOME")           # New raw data version
            i_start = lsearch_re(ls, "GENOME_START")
            i_end = lsearch_re(ls, "GENOME_STOP")
            header = "#Chrom\tstart\tend\tCOSMIC_ID\tCOSMIC_MUT_TYP"
            continue
        if line == "":
            continue
        cnv_id, cnv_type = li(ls, i_id), li(ls, i_type)
        if i_coord != -1:
            m = re.search(r"^([0-9]+):([0-9]+)\.\.([0-9]+)$", li(ls, i_coord))
            # Tcl parity: when the regexp fails, chrom/start/end keep the values of the previous line
            if m:
                chrom, start, end = m.group(1), m.group(2), m.group(3)
            if chrom is None:
                continue
        else:
            chrom, start, end = li(ls, i_chrom), li(ls, i_start), li(ls, i_end)
        to_write.setdefault(chrom, []).append(f"{chrom}\t{start}\t{end}\t{cnv_id}\t{cnv_type}")

    print(f"\t\t...creation of the COSMIC bed file in {users_dir}/")
    _write_text(header if header is not None else "", formatted)
    for ch in KARYOTYPE:
        if ch not in to_write:
            continue
        rows = sorted(set(to_write[ch]))
        rows = tcl.sort_on_element(tcl.sort_on_element(rows, 2), 1)
        _write_text("\n".join(rows), formatted)
    _delete(downloaded)
    print(f"\t\t   done ({_now()})")
    c["COSMICann"] = "1"


# --------------------------------------------------------------------------
# User BED files
# --------------------------------------------------------------------------
def users_dir(c) -> str:
    return f"{_ann_dir(c)}/Users/{c['genomeBuild']}"


def check_users_bed(config) -> None:
    """checkUsersBED."""
    c = config
    n_cols = number_of_annotation_col(c)
    for user_bed in glob.glob(f"{users_dir(c)}/*/*.bed"):
        not_formatted = re.sub(r"(.formatted.sorted)?.bed$", ".bed", user_bed, count=1, flags=re.I)
        formatted = re.sub(r"(.formatted.sorted)?.bed$", ".formatted.bed", user_bed, count=1, flags=re.I)
        formatted_sorted = re.sub(r"(.formatted.sorted)?.bed$", ".formatted.sorted.bed", user_bed, count=1,
                                  flags=re.I)
        header_file = re.sub(r".bed$", ".header.tsv", not_formatted, count=1, flags=re.I)
        if os.path.exists(not_formatted):
            check_bed(not_formatted, c["outputDir"], os.path.dirname(not_formatted))
            _unix_sort(c, formatted, formatted_sorted, "checkUsersBED")
            _delete(formatted)

        all_names = first_line_from_file(header_file).split("\t")
        ann_names = all_names[3:]
        n_header = len(all_names)
        n_bed = len(first_line_from_file(formatted_sorted).split("\t"))
        if n_header != n_bed or n_header <= 3:
            if n_header != n_bed:
                print(f"WARNING: {os.path.basename(header_file)} has {n_header} fields, but {n_bed} were expected. "
                      "Header bed file not used.\n")
            else:
                print(f"WARNING: {os.path.basename(header_file)} has {n_header} fields. More than 3 are expected. "
                      "Header bed file not used.\n")
            os.replace(header_file, header_file + ".bed")
            _write_text("\t".join("" for _ in range(n_bed)), header_file)
            # Tcl parity: L_ann_ColName_from_header becomes the number $nColBed-3, a one-element list
            ann_names = [str(n_bed - 3)]
        n_cols[formatted_sorted] = len(ann_names)


class UserBED:
    """userBEDannotation (also used for the cytoband)."""

    def __init__(self, config):
        self.config = config
        self.text: Dict[str, Dict[str, str]] = {}
        self.empty: Dict[str, str] = {}

    def _load(self, bed: str, use3points: bool) -> None:
        c = self.config
        n = number_of_annotation_col(c)[bed]
        self.empty[bed] = "\t".join("" for _ in range(n))
        tmp = _out_tmp(c, r".formatted.sorted.bed$", ".intersect.userBED")
        _delete(tmp)
        full_split = c["fullAndSplitBedFile"]
        bt = c.get("bedtools", "bedtools")
        prox = int(tcl.to_number(str(c.get("breakpointProximity", "100"))))
        if re.search("FtIncludedInSV|SVincludedInFt|AnyOverlap", bed):
            _intersect(c, full_split, bed, tmp, "userBEDannotation")
        elif "BNDproximity" in bed:
            ok, err = _exec_ok([bt, "window", "-a", full_split, "-b", bed, "-w", str(c["breakpointProximity"])], tmp)
            if not ok:
                raise ConfigExit(f"-- userBEDannotation --\n{bt} window -a {full_split} -b {bed} "
                                 f"-w {c['breakpointProximity']} > {tmp}\n{err}\nExit with error", 2)
        overlap = float(tcl.to_number(str(c.get("overlap", "100"))))
        reciprocal = str(c.get("reciprocal", "0")).lower() in ("1", "yes", "true", "on")
        ann: Dict[str, List[List[str]]] = {}
        if os.path.exists(tmp):
            lines = _read_lines(tmp)
        else:
            lines = []
        for line in lines:
            if line == "":
                continue
            ls = line.split("\t")
            sv_start, sv_end = int(ls[1]), int(ls[2])
            ft_end = int(ls[len(ls) - 1 - n])
            ft_start = int(ls[len(ls) - 2 - n])
            ft_length = ft_end - ft_start if ft_end - ft_start >= 1 else 1
            sv_length = sv_end - sv_start if sv_end - sv_start >= 1 else 1
            overlap_start = ft_start if sv_start < ft_start else sv_start
            overlap_end = sv_end if sv_end < ft_end else ft_end
            overlap_length = overlap_end - overlap_start
            if "FtIncludedInSV" in bed:
                if overlap_length * 100.0 / ft_length < overlap:
                    continue
            elif "SVincludedInFt" in bed:
                if overlap_length * 100.0 / sv_length < overlap:
                    continue
                if reciprocal and overlap_length * 100.0 / ft_length < overlap:
                    continue
            elif "AnyOverlap" in bed:
                pass
            elif "BNDproximity" in bed:
                if abs(sv_start - ft_start) > prox:
                    continue
                if abs(sv_end - ft_end) > prox:
                    continue
            key = f"{ls[0]},{ls[1]},{ls[2]}"
            cols = ann.setdefault(key, [[] for _ in range(n)])
            for k in range(n):
                cols[k].append(ls[len(ls) - n + k])
        texts: Dict[str, str] = {}
        for key, cols in ann.items():
            parts = []
            for values in cols:
                values = remove_redundancy_without_sorting(values)
                if use3points and len(values) > 50:
                    values = values[0:51] + ["..."]
                parts.append(";".join(values))
            texts[key] = "\t".join(parts)
        self.text[bed] = texts
        _delete(tmp)

    def __call__(self, bed: str, chrom: str, start: str, end: str, use3points: bool = True) -> str:
        if bed not in self.text:
            self._load(bed, use3points)
        return self.text[bed].get(f"{chrom},{start},{end}", self.empty[bed])


# --------------------------------------------------------------------------
# OrganizeAnnotation
# --------------------------------------------------------------------------
USER_KINDS = ("SVincludedInFt", "FtIncludedInSV", "AnyOverlap", "BNDproximity")


class RegionAnnotations:
    """The region columns of OrganizeAnnotation.

    ``cytoband_header`` goes right after ``Annotation_mode``. ``header_parts``
    returns the other parts, keyed by where they go:
    ``SVincludedInFt``, ``TAD``, ``FtIncludedInSV``, ``AnyOverlap``,
    ``BNDproximity`` and ``breakpoints`` (GC content, Repeat, Gap, SegDup,
    ENCODE blacklist).
    """

    def __init__(self, config):
        c = config
        self.config = c
        self.user_bed = UserBED(c)
        self.tad = TADannotation(c)
        self.gc = GCcontent(c)
        self.bp = {name: BreakpointAnnotation(c, name) for name in ("Repeat", "Gap", "SegDup", "ENCODEblacklist")}
        self.user_files: Dict[str, List[str]] = {}

    def _user_files(self, kind: str) -> List[str]:
        return glob.glob(f"{users_dir(self.config)}/{kind}/*.formatted.sorted.bed")

    def cytoband_header(self) -> str:
        c = self.config
        if _flag(c, "cytoband"):
            if "CytoBand" not in _wanted(c):
                c["cytoband"] = "0"
            else:
                return "\tCytoBand"
        return ""

    def header_parts(self) -> Dict[str, str]:
        c = self.config
        parts: Dict[str, str] = {}
        for kind in USER_KINDS:
            files = self._user_files(kind)
            self.user_files[kind] = files
            names = []
            for bed in files:
                header_file = re.sub(r".formatted.sorted.bed$", ".header.tsv", bed, count=1, flags=re.I)
                names.append("\t".join(first_line_from_file(header_file).split("\t")[3:]))
            parts[kind] = "\t" + "\t".join(names) if files else ""
        parts["TAD"] = self.tad.header()
        bp = ""
        if _flag(c, "gcContentAnn"):
            bp += "\tGC_content_left\tGC_content_right"
        for name in ("Repeat", "Gap", "SegDup", "ENCODEblacklist"):
            flag, cols, _, _ = BREAKPOINT_SOURCES[name]
            if _flag(c, flag):
                bp += "".join(f"\t{col}" for col in cols)
        parts["breakpoints"] = bp
        return parts

    def cytoband_text(self, chrom: str, start: str, end: str) -> str:
        c = self.config
        if not _flag(c, "cytoband"):
            return None
        bed = f"{cytoband_dir(c)}/cytoBand_{c['genomeBuild']}.formatted.sorted.bed"
        if not os.path.exists(bed):
            return ""
        bands = tcl.split_chars(self.user_bed(bed, chrom, start, end, use3points=False), ";")
        if not bands:
            return ""
        text = bands[0]
        if len(bands) > 1:
            text += f"-{bands[-1]}"
        return text

    def texts(self, chrom: str, start: str, end: str, sv_left: str, sv_right: str, mode: str) -> Dict[str, str]:
        """``start``/``end``: the SV (full lines) or the intersection (split lines)."""
        c = self.config
        out: Dict[str, str] = {}
        for kind in USER_KINDS:
            files = self.user_files.get(kind, [])
            out[kind] = "\t" + "\t".join(self.user_bed(bed, chrom, start, end) for bed in files) if files else ""
        if _flag(c, "tadAnn"):
            out["TAD"] = "\t" + self.tad(chrom, start, end)
        else:
            out["TAD"] = ""
        bp = ""
        full = mode == "full"
        if _flag(c, "gcContentAnn"):
            bp += "\t" + (f"{self.gc(chrom, sv_left)}\t{self.gc(chrom, sv_right)}" if full else "\t")
        for name in ("Repeat", "Gap", "SegDup", "ENCODEblacklist"):
            flag = BREAKPOINT_SOURCES[name][0]
            if _flag(c, flag):
                ann = self.bp[name]
                if full:
                    bp += f"\t{ann(chrom, sv_left)}\t{ann(chrom, sv_right)}"
                else:
                    bp += "\t" + ("\t\t\t" if ann.two_cols else "\t")
        out["breakpoints"] = bp
        return out
