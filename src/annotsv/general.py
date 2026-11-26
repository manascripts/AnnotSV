"""Port of ``share/tcl/AnnotSV/AnnotSV-general.tcl``.

File helpers, AnnotSV_ID building, SV type normalisation, annotation merging
and BED file checks. Function names follow PEP 8; the Tcl proc each one
replaces is named in its docstring.
"""

from __future__ import annotations

import gzip
import os
import re
from typing import Dict, Iterable, List, Optional, Tuple

from . import tcl

# --------------------------------------------------------------------------
# Sort comparators (AscendingSortOnElementN / DescendingSortOnElementN)
# --------------------------------------------------------------------------
sort_on_element = tcl.sort_on_element


# --------------------------------------------------------------------------
# Working with files
# --------------------------------------------------------------------------
def _open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return open(path, "r", newline="")


def first_line_from_file(path: str) -> str:
    """FirstLineFromFile: first line of a plain or ``.gz`` file."""
    with _open_text(path) as f:
        line = f.readline()
    return line.rstrip("\n").rstrip("\r") if line.endswith("\r\n") else line.rstrip("\n")


def content_from_file(path: str) -> str:
    """ContentFromFile: whole file content, without the final newline."""
    if path == "":
        return ""
    with open(path, "r", newline="") as f:
        text = f.read()
    return text[:-1] if text.endswith("\n") else text


def lines_from_file(path: str) -> List[str]:
    """LinesFromFile: ``split [ContentFromFile] "\\n"``."""
    return content_from_file(path).split("\n")


def content_from_gz_file(path: str) -> str:
    """ContentFromGZFile."""
    if path == "":
        return ""
    with gzip.open(path, "rt", newline="") as f:
        text = f.read()
    return text[:-1] if text.endswith("\n") else text


def lines_from_gz_file(path: str) -> List[str]:
    """LinesFromGZFile."""
    return content_from_gz_file(path).split("\n")


def replace_text_in_file(text: str, path: str) -> int:
    """ReplaceTextInFile: overwrite ``path`` with ``text`` plus a newline."""
    with open(path, "w", newline="") as f:
        f.write(text + "\n")
    return 1


def write_text_in_file(text: str, path: str) -> int:
    """WriteTextInFile: append ``text`` plus a newline to ``path``."""
    with open(path, "a", newline="") as f:
        f.write(text + "\n")
    return 1


# --------------------------------------------------------------------------
# AnnotSV_ID
# --------------------------------------------------------------------------
_BRACKET_ALT = re.compile(r"([acgtnACGTN]*)(\[|\])([^:]+):([0-9]+)(\[|\])([acgtnACGTN]*)")


def remove_ref_in_bracketed_alt(ref: str, alt: str) -> Tuple[str, str]:
    """removeREFinBracketedALT.

    For square-bracketed ALT only: drop the REF base(s) from the ALT, replace
    every ``[`` / ``]`` with ``;`` and blank the REF. Returns ``(ref, alt)``.
    """
    m = _BRACKET_ALT.search(alt)
    if not m:
        return ref, alt
    base_left, bracket_left, chrom, start, bracket_right, base_right = m.groups()
    if len(base_left) > 1:
        base_left = re.sub("^" + ref, "", base_left, count=1, flags=re.I)
    elif len(base_right) > 1:
        base_right = re.sub(ref + "$", "", base_right, count=1, flags=re.I)
    else:
        base_left = base_right = ""
    alt = f"{base_left}{bracket_left}{chrom}:{start}{bracket_right}{base_right}"
    alt = re.sub(r"\]|\[", ";", alt)
    return "", alt


class AnnotSVIDRegistry:
    """settingOfTheAnnotSVID, with the ``g_ID`` / ``g_deb`` globals as state."""

    def __init__(self) -> None:
        self._ids: Dict[Tuple[str, str, str], str] = {}
        self._used: set = set()

    def get(self, deb: str, ref: str, alt: str) -> str:
        """``deb`` is ``"${chrom}_${pos}_${end}_${svtype}"``."""
        ref, alt = remove_ref_in_bracketed_alt(ref, alt)
        key = (deb, ref, alt)
        if key not in self._ids:
            i = 1
            while f"{deb}_{i}" in self._used:
                i += 1
            self._used.add(f"{deb}_{i}")
            self._ids[key] = f"{deb}_{i}"
        return self._ids[key]


# --------------------------------------------------------------------------
# rsID
# --------------------------------------------------------------------------
def is_not_an_rs(rs: str) -> int:
    """isNotAnRS: 0 if ``rs`` looks like ``rs123``, else 1."""
    return 0 if re.search(r"^rs[0-9]+$", rs) else 1


# --------------------------------------------------------------------------
# Annotation lists
# --------------------------------------------------------------------------
def merge_annotation(l_ann: List[str], number_of_ann_col: int) -> str:
    """MergeAnnotation.

    ``l_ann`` holds one tab-separated annotation string per gene; the result
    joins each column across genes with ``;`` (empty values collapsed).
    """
    if len(l_ann) == 1:
        # Tcl: [join $L_ann] on a one-element list returns the element itself
        return l_ann[0]
    table: Dict[Tuple[int, int], str] = {}
    for g, gene_ann in enumerate(l_ann):
        if gene_ann == "":
            table[(g, 0)] = ""
        else:
            for c, col_ann in enumerate(gene_ann.split("\t")):
                table[(g, c)] = col_ann
    new_ann = []
    for c in range(number_of_ann_col):
        # Tcl would raise on a missing (g, c) cell; we treat it as empty.
        one = ";".join(table.get((g, c), "") for g in range(len(l_ann)))
        one = re.sub(";;+", ";", one)
        one = re.sub("^;", "", one, count=1)
        one = re.sub(";$", "", one, count=1)
        new_ann.append(one)
    return "\t".join(new_ann)


def remove_redundancy_without_sorting(items: Iterable[str]) -> List[str]:
    """RemoveRedundancyWithoutSorting."""
    seen: List[str] = []
    for item in items:
        if item not in seen:
            seen.append(item)
    return seen


# --------------------------------------------------------------------------
# SV type
# --------------------------------------------------------------------------
def normalize_svtype(svtype: str) -> str:
    """normalizeSVtype: map a raw SV type to DEL, DUP, INV, INS, TRA or ""."""
    if re.search(r"del|loss|<CN0>|<CN1>", svtype, re.I):
        return "DEL"
    if re.search(r"dup|gain|MCNV", svtype, re.I):
        return "DUP"
    m = re.search(r"<CN([0-9]+)>", svtype, re.I)
    if m:
        # Tcl keeps the original value when <CNn> has n <= 1 (unreachable for
        # CN0/CN1, already caught above).
        return "DUP" if int(m.group(1)) > 1 else svtype
    if re.search(r"inv", svtype, re.I):
        return "INV"
    if re.search(r"ins|MEI|alu|line|sva", svtype, re.I):
        return "INS"
    if re.search(r"TRA|TRN", svtype, re.I):
        return "TRA"
    return ""


# --------------------------------------------------------------------------
# BED files
# --------------------------------------------------------------------------
_KARYOTYPE = [str(i) for i in range(1, 23)] + ["X", "Y", "M", "MT"]


def check_bed(bed_file: str, output_dir: str, bed_dir: str = "", merge_overlap: bool = False):
    """checkBed.

    Writes ``*.formatted.bed`` (sorted, "chr" removed, duplicates removed,
    optionally merged) into ``bed_dir`` or ``output_dir``, then deletes the
    input file. Returns the new path, or ``1`` if it already existed (as in
    Tcl). For BED files under a ``/Users/`` directory a ``*.header.tsv`` is
    created next to the input when missing.
    """
    base = re.sub(r"\.bed$", ".formatted.bed", os.path.basename(bed_file), flags=re.I)
    new_bed = os.path.join(bed_dir if bed_dir else output_dir, base)
    if os.path.exists(new_bed):
        return 1

    text: Dict[str, List[str]] = {}
    lengths: List[int] = []
    header = ""
    with open(bed_file, "r", newline="") as f:
        for raw in f:
            line = raw.rstrip("\n")
            if line == "":
                continue
            line = "\t".join(v.strip() for v in line.split("\t"))
            if line.startswith("#"):
                header = line
                continue
            first = tcl.lindex(line, 0)
            if "chr" in first:
                chrom = first.replace("chr", "", 1)
                text.setdefault(chrom, []).append(line[3:])
            else:
                text.setdefault(first, []).append(line)
            lengths.append(len(line.split("\t")))

    unique_lengths = sorted(set(lengths))
    if len(unique_lengths) != 1:
        print(f"{bed_file}: Unexpected file format, not the same length for each line.")
        print("\t-> not using the associated annotations")
        use_ann = False
    else:
        use_ann = unique_lengths[0] > 3
    # Tcl keeps a list here when lengths differ (and then fails); use the max.
    n_cols = max(unique_lengths) if unique_lengths else 0

    if "/Users/" in bed_file:
        user_header = re.sub(r"(\.formatted)?\.bed$", ".header.tsv", bed_file, flags=re.I)
        if not os.path.exists(user_header):
            if len(header.split("\t")) == n_cols:
                write_text_in_file(header, user_header)
            else:
                write_text_in_file("\t".join([""] * n_cols), user_header)

    for chrom in _KARYOTYPE:
        if chrom not in text:
            continue
        lines = sorted(set(text[chrom]))
        lines = sort_on_element(sort_on_element(lines, 2), 1)

        if merge_overlap:
            new_text = []
            ls = lines[0].split("\t")
            prev_start, prev_end = ls[1], ls[2]
            prev_ann = {j: _get(ls, j) for j in range(3, n_cols)} if use_ann else {}
            ann: Dict[int, str] = {}
            for line in lines[1:]:
                ls = line.split("\t")
                start, end = ls[1], ls[2]
                if use_ann:
                    ann = {j: _get(ls, j) for j in range(3, n_cols)}
                if tcl.compare(start, prev_end) < 0:
                    if tcl.compare(end, prev_end) > 0:
                        prev_end = end
                    if use_ann:
                        for j in range(3, n_cols):
                            prev_ann[j] = f"{prev_ann[j]}/{ann[j]}"
                else:
                    new_text.append(_merged_line(chrom, prev_start, prev_end, prev_ann, use_ann, n_cols))
                    prev_start, prev_end = start, end
                    if use_ann:
                        prev_ann = dict(ann)
            new_text.append(_merged_line(chrom, prev_start, prev_end, prev_ann, use_ann, n_cols))
            write_text_in_file("\n".join(new_text), new_bed)
        else:
            write_text_in_file("\n".join(lines), new_bed)
        del text[chrom]

    os.remove(bed_file)
    return new_bed


def _get(ls: List[str], j: int) -> str:
    return ls[j] if j < len(ls) else ""


def _merged_line(chrom, start, end, ann, use_ann, n_cols) -> str:
    out = f"{chrom}\t{start}\t{end}"
    if use_ann:
        out += "".join(f"\t{ann[j]}" for j in range(3, n_cols))
    return out


def is_an_empty_file(path: str) -> int:
    """isAnEmptyFile: 1 if not a .bed/.vcf(.gz) or no non-comment line."""
    if not re.search(r"\.(bed|vcf(.gz)?)$", path, re.I):
        return 1
    opener = gzip.open if re.search(r".gz$", path, re.I) else open
    with opener(path, "rt") as f:
        for line in f:
            line = line.rstrip("\n")
            if line == "" or line.startswith("#"):
                continue
            return 0
    return 1


def _first_data_line_length(bed_file: str) -> Optional[int]:
    with open(bed_file, "r", newline="") as f:
        for raw in f:
            line = raw.rstrip("\n")
            if line == "" or line.startswith("#"):
                continue
            return len(line.split("\t"))
    return None


def check_samplesid_bed_col(config, bed_file: str) -> None:
    """checksamplesidBEDcol: validate ``-samplesidBEDcol`` and normalise the
    sample separators in that column to commas (file rewritten in place)."""
    if int(config["samplesidBEDcol"]) == -1:
        return
    length = _first_data_line_length(bed_file) or 0
    col = int(config["samplesidBEDcol"])
    if col > length:
        print(f"\nWARNING: -samplesidBEDcol = {config['samplesidBEDcol']}")
        print(f"This value should correspond to a column of {bed_file} ({length} columns)")
        print("-samplesidBEDcol is set to -1.\n")
        config["samplesidBEDcol"] = "-1"
    elif col < 4:
        print(f"\nWARNING: -samplesidBEDcol = {config['samplesidBEDcol']}")
        print("This value should be > 3")
        print("-samplesidBEDcol is set to -1.\n")
        config["samplesidBEDcol"] = "-1"

    if int(config["samplesidBEDcol"]) == -1:
        return
    idx = int(config["samplesidBEDcol"]) - 1
    config["samplesidBEDcol"] = str(idx)
    config["samplesidTSVcol"] = str(idx + 2)

    to_write, from_start = [], []
    with open(bed_file, "r", newline="") as f:
        for raw in f:
            line = raw.rstrip("\n")
            if line == "":
                continue
            from_start.append(line)
            if line.startswith("#"):
                to_write.append(line)
                continue
            ls = line.split("\t")
            part2 = _get(ls, idx)
            # Tcl: [split $part2 " +| *\| *| *, *| *; *| */ *"] splits on each
            # character of that string: space + | * , ; /
            part2 = ",".join(tcl.split_chars(part2, " +|*,;/"))
            if idx < len(ls):
                ls[idx] = part2
            else:
                ls.append(part2)
            to_write.append("\t".join(ls))
    if to_write != from_start:
        replace_text_in_file("\n".join(to_write), bed_file)


def check_svt_bed_col(config, bed_file: str) -> None:
    """checksvtBEDcol: validate ``-svtBEDcol`` (disables ranking if invalid)."""
    if int(config["svtBEDcol"]) == -1:
        return
    length = _first_data_line_length(bed_file) or 0
    col = int(config["svtBEDcol"])
    if col > length:
        print(f"\nWARNING: -svtBEDcol = {config['svtBEDcol']}")
        print(f"This value should correspond to a column of {bed_file} ({length} columns)")
        print("-svtBEDcol is set to -1. NO SV RANKING WILL BE PERFORMED.\n")
        config["svtBEDcol"] = "-1"
        config["ranking"] = "0"
    elif col < 4:
        print(f"\nWARNING: -svtBEDcol = {config['svtBEDcol']}")
        print("This value should be > 3")
        print("-svtBEDcol is set to -1. NO SV RANKING WILL BE PERFORMED.\n")
        config["svtBEDcol"] = "-1"
        config["ranking"] = "0"
    else:
        config["svtBEDcol"] = str(col - 1)
        config["svtTSVcol"] = str(col + 1)
