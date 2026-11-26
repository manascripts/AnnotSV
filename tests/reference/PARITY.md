# Output parity policy: Python port vs. Tcl AnnotSV

**Standard: exact text, with named exceptions.** A Python run passes when
`compare_outputs.py GOLDEN CANDIDATE` exits 0 with its default options.

## What must match exactly

* **Output files**: the same set of files, at the same relative paths.
* **TSV header**: the same columns, in the same order.
* **Rows**: the same rows, in the same order. Tcl sorts by ranking score, then
  by Exomiser score, then by BED order (`AnnotSV-write.tcl`, a stable
  `lsort -command`). Python's `sorted(..., reverse=True)` is also stable and
  gives the same order when it gets the same input order.
* **Cell text**: compared as strings, with no numeric tolerance. Most floats are
  written with an explicit `format "%.Nf"` (use `f"{x:.Nf}"`). Unformatted
  doubles print the shortest way in Tcl 8.5+, the same as Python's `repr`.
  Integers must stay integers: print `5`, not `5.0`.
* **VCF and other text outputs**: compared line for line.

## Tcl quirks are kept

Reproduce Tcl's behavior even when it looks like a bug, and mark each case with
a `# Tcl parity:` comment. Fixing one is a separate, deliberate change that also
updates the golden files. Known cases:

* **Integer division.** In `expr {$a*100/$b}`, `/` floors when both operands are
  integers (for example, `AnnotSV-genes.tcl` `percentCDS` and the
  `AnnotSV-benignsv.tcl` frequencies). In Python, use `//` with Tcl's floor
  semantics. Tcl and Python both floor toward −∞.
* **String sort of numbers.** A plain `lsort` sorts ASCII strings, so `"10"` comes
  before `"9"` (for example, the min LOEUF/pLI in `AnnotSV-loeuf-pLI.tcl`).
  Use `sorted()` on the strings, not on floats.
* **`expr` comparisons** inside the `*SortOnElement*` procs compare numerically
  when both values look numeric, and as strings otherwise.

## Named exceptions

These are the only relaxations. They are on by default in `compare_outputs.py`,
and `--strict` turns them off. An entry belongs here only when the Tcl output
order comes from hash iteration (`array names`) and nothing sorts it afterwards.
Each entry must name its source.

| What | Relaxation | Tcl source |
|---|---|---|
| `*.unannotated.tsv` | lines compared as a set | `array names TextToWrite_rescue` in `VCFsToBED` (AnnotSV-vcf.tcl) writes the reciprocal-breakend notes |
| list columns | none so far | The `array names` loops audited (OMIM, regulatory elements, TAD, user BED, ExAC, benign SV) all `lsort -unique` their lists or keep the input-file order before writing them |

If a new mismatch shows up, trace it to an `array names` loop before you add a
column with `DEFAULT_UNORDERED_COLS` or a file with
`DEFAULT_UNORDERED_LINE_FILES`. Otherwise, fix the port.

## Not compared

* `*.log` files (timestamps, absolute paths). Use `--include-logs` to diff them anyway.
* Annotation cache files that AnnotSV builds under `share/AnnotSV/Annotations_*`
  (for example, the ExAC CNV-intolerance and benign-SV gene files, which are
  written in hash order and only used for lookups).
* Exomiser scores also depend on the Java/Exomiser install. If they differ, check
  the Exomiser version first.

`--rel-tol`, `--abs-tol` and `--unordered-lists` remain available for
debugging. A run that passes only with them does not meet the standard.
