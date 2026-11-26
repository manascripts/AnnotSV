#!/usr/bin/env python3
"""Column-by-column comparison of AnnotSV outputs against the Tcl golden reference.

Compares either two files or two directory trees (files are paired by relative path).

  python3 tests/reference/compare_outputs.py GOLDEN CANDIDATE [options]

  # whole golden tree vs. a Python-port run laid out the same way
  python3 tests/reference/compare_outputs.py tests/reference/golden /tmp/py_run

  # one file, ignoring list order inside ';'-separated cells
  python3 tests/reference/compare_outputs.py a.tsv b.tsv --unordered-lists

*.tsv files are compared as tables: header (column set and order), row set
(rows are matched by key columns, not by position), then every shared column
cell by cell. Other text files (VCF, logs, ...) are compared line by line;
*.log files are skipped unless --include-logs is given.

Exit status: 0 = identical (under the chosen tolerances), 1 = differences, 2 = usage error.
"""
from __future__ import annotations

import argparse
import difflib
import gzip
import math
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

# Columns that together identify an output row of an AnnotSV .annotated.tsv.
# AnnotSV_ID alone is shared by the "full" row and its "split" rows.
DEFAULT_KEY = ["AnnotSV_ID", "Annotation_mode", "Gene_name", "Tx"]


# --------------------------------------------------------------------------- io
def read_text(path: Path) -> str:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace", newline="") as fh:
        return fh.read()


def read_tsv(path: Path) -> tuple[list[str], list[list[str]]]:
    lines = read_text(path).split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    lines = [ln[:-1] if ln.endswith("\r") else ln for ln in lines]
    if not lines:
        return [], []
    header = lines[0].split("\t")
    rows = [ln.split("\t") for ln in lines[1:]]
    return header, rows


# ----------------------------------------------------------------- comparison
@dataclass
class Options:
    key: list[str]
    rel_tol: float
    abs_tol: float
    unordered_lists: bool
    list_seps: str
    ignore_cols: set[str]
    max_examples: int
    include_logs: bool


def _as_float(s: str):
    try:
        f = float(s)
    except ValueError:
        return None
    return None if math.isnan(f) else f


def _norm_list(s: str, seps: str) -> str:
    for sep in seps:
        if sep in s:
            return sep.join(sorted(s.split(sep)))
    return s


def cells_equal(a: str, b: str, opt: Options) -> bool:
    if a == b:
        return True
    if opt.rel_tol or opt.abs_tol:
        fa, fb = _as_float(a), _as_float(b)
        if fa is not None and fb is not None and math.isclose(fa, fb, rel_tol=opt.rel_tol, abs_tol=opt.abs_tol):
            return True
    if opt.unordered_lists and _norm_list(a, opt.list_seps) == _norm_list(b, opt.list_seps):
        return True
    return False


@dataclass
class FileReport:
    name: str
    lines: list[str] = field(default_factory=list)
    ok: bool = True

    def fail(self, msg: str) -> None:
        self.ok = False
        self.lines.append(msg)

    def info(self, msg: str) -> None:
        self.lines.append(msg)


def _short(s: str, n: int = 80) -> str:
    return s if len(s) <= n else s[: n - 3] + "..."


def compare_tsv(gold: Path, cand: Path, opt: Options, rep: FileReport) -> None:
    gh, grows = read_tsv(gold)
    ch, crows = read_tsv(cand)

    # 1. header
    if gh != ch:
        missing = [c for c in gh if c not in ch]
        extra = [c for c in ch if c not in gh]
        if missing:
            rep.fail(f"columns missing from candidate ({len(missing)}): {', '.join(missing)}")
        if extra:
            rep.fail(f"extra columns in candidate ({len(extra)}): {', '.join(extra)}")
        common_g = [c for c in gh if c in ch]
        common_c = [c for c in ch if c in gh]
        if common_g != common_c:
            rep.fail("column order differs: first mismatch at "
                     + next(f"#{i + 1} gold={a!r} cand={b!r}" for i, (a, b) in enumerate(zip(common_g, common_c)) if a != b))
        dup = [c for c, n in Counter(gh).items() if n > 1]
        if dup:
            rep.info(f"note: duplicate column names in golden header: {dup}")

    gidx = {c: i for i, c in enumerate(gh)}
    cidx = {c: i for i, c in enumerate(ch)}
    shared = [c for c in gh if c in cidx and c not in opt.ignore_cols]

    # 2. rows: key by the key columns present in both files, else by position
    key_cols = [k for k in opt.key if k in gidx and k in cidx]

    def keyed(rows, idx):
        out: dict[tuple, list[str]] = {}
        seen: Counter = Counter()
        for r in rows:
            base = tuple(r[idx[k]] if idx[k] < len(r) else "" for k in key_cols) if key_cols else ()
            seen[base] += 1
            out[base + (seen[base],)] = r  # occurrence counter keeps duplicates distinct
        return out

    if key_cols:
        gk, ck = keyed(grows, gidx), keyed(crows, cidx)
    else:
        gk = {(i,): r for i, r in enumerate(grows)}
        ck = {(i,): r for i, r in enumerate(crows)}
        rep.info("note: no key columns found, rows matched by position")

    if len(grows) != len(crows):
        rep.fail(f"row count differs: gold={len(grows)} cand={len(crows)}")
    only_g = [k for k in gk if k not in ck]
    only_c = [k for k in ck if k not in gk]
    if only_g:
        rep.fail(f"{len(only_g)} row(s) only in golden, e.g. " + "; ".join(map(_fmt_key, only_g[: opt.max_examples])))
    if only_c:
        rep.fail(f"{len(only_c)} row(s) only in candidate, e.g. " + "; ".join(map(_fmt_key, only_c[: opt.max_examples])))

    common = [k for k in gk if k in ck]
    if key_cols:
        gorder = [k for k in gk if k in ck]
        corder = [k for k in ck if k in gk]
        if gorder != corder:
            rep.fail("row order differs (same rows, different order)")

    # 3. cells, column by column
    diffs: dict[str, list[tuple]] = defaultdict(list)
    for k in common:
        gr, cr = gk[k], ck[k]
        for col in shared:
            gi, ci = gidx[col], cidx[col]
            a = gr[gi] if gi < len(gr) else ""
            b = cr[ci] if ci < len(cr) else ""
            if not cells_equal(a, b, opt):
                diffs[col].append((k, a, b))
    for col in shared:
        if col in diffs:
            d = diffs[col]
            rep.fail(f"column {col!r}: {len(d)}/{len(common)} cell(s) differ")
            for k, a, b in d[: opt.max_examples]:
                rep.info(f"    {_fmt_key(k)}: gold={_short(a)!r} cand={_short(b)!r}")
    if rep.ok:
        rep.info(f"identical: {len(gh)} columns x {len(grows)} rows")


def _fmt_key(k: tuple) -> str:
    *parts, n = k if len(k) > 1 else (*k, 1)
    s = "|".join(str(p) for p in parts) or "<row>"
    return s + (f" (#{n})" if n != 1 else "")


def compare_text(gold: Path, cand: Path, opt: Options, rep: FileReport) -> None:
    g = read_text(gold).splitlines()
    c = read_text(cand).splitlines()
    if g == c:
        rep.info(f"identical: {len(g)} lines")
        return
    diff = list(difflib.unified_diff(g, c, "golden", "candidate", lineterm="", n=0))
    rep.fail(f"text differs ({sum(1 for l in diff if l[:1] in '+-' and l[:3] not in ('+++', '---'))} changed lines)")
    for line in diff[2: 2 + 4 * opt.max_examples]:
        rep.info("    " + _short(line, 160))


def compare_file(gold: Path, cand: Path, opt: Options, name: str) -> FileReport:
    rep = FileReport(name)
    if not cand.exists():
        rep.fail("missing in candidate")
        return rep
    stem = gold.name[:-3] if gold.name.endswith(".gz") else gold.name
    try:
        if stem.endswith(".tsv"):
            compare_tsv(gold, cand, opt, rep)
        else:
            compare_text(gold, cand, opt, rep)
    except Exception as exc:  # keep going over the other files
        rep.fail(f"could not compare: {exc!r}")
    return rep


SKIP_NAMES = {"SUMMARY.tsv"}


def iter_files(root: Path, opt: Options):
    for p in sorted(root.rglob("*")):
        if not p.is_file() or p.name in SKIP_NAMES:
            continue
        if p.suffix == ".log" and not opt.include_logs:
            continue
        yield p.relative_to(root)


# ------------------------------------------------------------------------ cli
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("golden", type=Path, help="golden file or directory (Tcl reference)")
    ap.add_argument("candidate", type=Path, help="candidate file or directory (e.g. Python port)")
    ap.add_argument("--key", default=",".join(DEFAULT_KEY),
                    help="comma-separated row key columns (default: %(default)s); missing ones are dropped")
    ap.add_argument("--rel-tol", type=float, default=0.0, help="relative tolerance for numeric cells")
    ap.add_argument("--abs-tol", type=float, default=0.0, help="absolute tolerance for numeric cells")
    ap.add_argument("--unordered-lists", action="store_true",
                    help="treat cells holding ';' or ',' separated lists as unordered")
    ap.add_argument("--list-seps", default=";,", help="separators used by --unordered-lists (default: ';,')")
    ap.add_argument("--ignore-col", action="append", default=[], help="column to skip (repeatable)")
    ap.add_argument("--max-examples", type=int, default=5, help="example cells shown per differing column")
    ap.add_argument("--include-logs", action="store_true", help="also diff *.log files")
    ap.add_argument("-q", "--quiet", action="store_true", help="only print files that differ, and the summary")
    a = ap.parse_args(argv)

    opt = Options(key=[k for k in a.key.split(",") if k], rel_tol=a.rel_tol, abs_tol=a.abs_tol,
                  unordered_lists=a.unordered_lists, list_seps=a.list_seps, ignore_cols=set(a.ignore_col),
                  max_examples=a.max_examples, include_logs=a.include_logs)

    if not a.golden.exists():
        ap.error(f"{a.golden} does not exist")
    if a.golden.is_file():
        cand = a.candidate / a.golden.name if a.candidate.is_dir() else a.candidate
        reports = [compare_file(a.golden, cand, opt, a.golden.name)]
        extra = []
    else:
        if not a.candidate.is_dir():
            ap.error(f"{a.candidate} is not a directory")
        rels = list(iter_files(a.golden, opt))
        reports = [compare_file(a.golden / r, a.candidate / r, opt, str(r)) for r in rels]
        extra = sorted(set(iter_files(a.candidate, opt)) - set(rels))

    n_bad = 0
    for rep in reports:
        if rep.ok and a.quiet:
            continue
        print(f"{'OK  ' if rep.ok else 'DIFF'} {rep.name}")
        for line in rep.lines:
            if not (rep.ok and line.startswith("identical")):
                print(f"     {line}")
        n_bad += not rep.ok
    for r in extra:
        print(f"NEW  {r}  (only in candidate)")

    print(f"\n{len(reports) - n_bad}/{len(reports)} file(s) match"
          + (f", {len(extra)} extra file(s) in candidate" if extra else ""))
    return 1 if n_bad else 0


if __name__ == "__main__":
    sys.exit(main())
