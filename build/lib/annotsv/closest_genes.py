"""Nearest genes: port of ``AnnotSV-closestgenes.tcl``.

``ClosestGenes(config)`` stands for ``closestGenesAnnotation``: called as
``closest(chrom, pos, "left"|"right")`` it returns the genes (in genes file
order, without duplicates) found up to 5 Mb to the left of ``pos`` (an SV
start) or to the right of ``pos`` (an SV end). The output writer keeps the
nearest one that the SV does not overlap.

As in Tcl, everything is computed on the first call, from the
``bedFile`` current at that time, with the same sort and ``bedtools
intersect`` steps.
"""

from __future__ import annotations

import os
import re
import subprocess
from typing import Dict, List

from .config import ConfigExit
from .general import write_text_in_file
from .genes import genes_dir, unix_sort_bed

WINDOW = 5000000
BATCH = 100000


def _tmp_name(config, suffix: str) -> str:
    c = config
    name = re.sub(r"(.formatted)?.bed$", lambda _m: suffix, c["bedFile"], count=1, flags=re.I)
    return f"{c['outputDir']}/{os.path.basename(name)}"


def _remove(path: str) -> None:
    if os.path.isfile(path):
        os.remove(path)


def _intersect(config, a_file: str, b_file: str, out_file: str, side: str) -> None:
    c = config
    base = [c.get("bedtools", "bedtools"), "intersect"]
    tail = ["-a", a_file, "-b", b_file, "-wa", "-wb"]
    for cmd in (base + ["-sorted"] + tail, base + tail):
        with open(out_file, "w") as out:
            r = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE)
        # Tcl "exec" fails on a non-zero exit status or on any stderr output
        if r.returncode == 0 and not r.stderr:
            return
    raise ConfigExit(f"-- closestGenesAnnotation, {side} intersect --\n"
                     f"{' '.join(base + ['-sorted'] + tail)} > {out_file}\n"
                     f"{r.stderr.decode(errors='replace')}\nExit with error", 2)


def _load(path: str, pos_col: int) -> Dict[str, List[str]]:
    genes: Dict[str, List[str]] = {}
    with open(path, "r", newline="") as f:
        for raw in f:
            line = raw[:-1] if raw.endswith("\n") else raw
            if line == "":
                continue
            ls = line.split("\t")
            key = f"{ls[0]},{ls[pos_col]}"
            gene = ls[7] if len(ls) > 7 else ""
            lst = genes.setdefault(key, [])
            if gene not in lst:
                lst.append(gene)
    return genes


class ClosestGenes:
    """closestGenesAnnotation (``g_closestGenesLeft`` / ``g_closestGenesRight``)."""

    def __init__(self, config):
        self.config = config
        self.left: Dict[str, List[str]] = {}
        self.right: Dict[str, List[str]] = {}
        self.done = False

    def _compute(self) -> None:
        c = self.config
        left_file = _tmp_name(c, ".closestLeft.bed")
        right_file = _tmp_name(c, ".closestRight.bed")
        _remove(f"{left_file}.tmp")
        _remove(f"{right_file}.tmp")

        # SV left coordinates: "start-5000000bp" "start" => closest left gene
        # SV right coordinates: "end" "end+5000000bp" => closest right gene
        i = 0
        to_left: List[str] = []
        to_right: List[str] = []
        with open(c["bedFile"], "r", newline="") as f:
            for raw in f:
                line = raw[:-1] if raw.endswith("\n") else raw
                if line == "":
                    continue
                ls = line.split("\t")
                chrom = ls[0]
                start = ls[1] if len(ls) > 1 else ""
                end = ls[2] if len(ls) > 2 else ""
                start_left = max(int(start) - WINDOW, 0)
                end_right = int(end) + WINDOW
                to_left.append(f"{chrom}\t{start_left}\t{start}")
                to_right.append(f"{chrom}\t{end}\t{end_right}")
                i += 1
                if i > BATCH:
                    write_text_in_file("\n".join(sorted(set(to_left))), f"{left_file}.tmp")
                    write_text_in_file("\n".join(sorted(set(to_right))), f"{right_file}.tmp")
                    to_left, to_right, i = [], [], 0
        write_text_in_file("\n".join(sorted(set(to_left))), f"{left_file}.tmp")
        write_text_in_file("\n".join(sorted(set(to_right))), f"{right_file}.tmp")

        for tmp_file, side in ((left_file, "left"), (right_file, "right")):
            _remove(tmp_file)
            unix_sort_bed(f"{tmp_file}.tmp", tmp_file, f"closestGenesAnnotation, {side} sort")
            _remove(f"{tmp_file}.tmp")

        genes_file = f"{genes_dir(c)}/genes.{c['tx']}.sorted.bed"
        inter_left = _tmp_name(c, ".intersect.closestGenesLeft.bed")
        inter_right = _tmp_name(c, ".intersect.closestGenesRight.bed")
        for path in (inter_left, inter_right):
            _remove(path)
        _intersect(c, left_file, genes_file, inter_left, "left")
        _intersect(c, right_file, genes_file, inter_right, "right")

        self.left = _load(inter_left, 2)
        _remove(inter_left)
        self.right = _load(inter_right, 1)
        _remove(inter_right)

        self.done = True
        _remove(left_file)
        _remove(right_file)

    def __call__(self, chrom: str, pos: str, side: str) -> List[str]:
        if not self.done:
            self._compute()
        if side == "left":
            return list(self.left.get(f"{chrom},{pos}", []))
        if side == "right":
            return list(self.right.get(f"{chrom},{pos}", []))
        return []
