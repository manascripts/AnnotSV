# Tcl reference outputs (golden files)

Golden outputs of the original Tcl AnnotSV, used to check the Python port column by column.

**Pass/fail standard:** exact text with a few named exceptions. See [PARITY.md](PARITY.md).

| File | Purpose |
|---|---|
| `run_tcl_reference.sh` | Runs every `tests/AnnotSV/test_*/command_Tcl.sh` on a scratch copy of the inputs and copies all produced files to `golden/<test>/`, plus `golden/SUMMARY.tsv` (status, exit code, runtime per test). |
| `compare_outputs.py` | Compares two files or two directory trees. `.tsv` files: header, row set (keyed on `AnnotSV_ID, Annotation_mode, Gene_name, Tx`), row order, then every column cell by cell. Other files: line diff. `*.log` skipped unless `--include-logs`. Exit 0 = match, 1 = differences. |
| `PARITY.md` | What must match, which Tcl quirks are kept, and the named exceptions. |
| `golden/` | The captured outputs (created by the run script). |

## Prerequisites

* `tclsh` + `tcllib`, `bedtools`, `bcftools`, `java` (Exomiser), `python3`
* Annotations installed into the repo: `make PREFIX=. install` then `make PREFIX=. install-human-annotation` and `make PREFIX=. install-mouse-annotation` (test_19). Needs network access to `www.lbgi.fr`, `data.monarchinitiative.org` and GitHub release downloads.

## Usage

```bash
export ANNOTSV=$PWD                                   # repo root
tests/reference/run_tcl_reference.sh                   # all tests (~1h+)
tests/reference/run_tcl_reference.sh 'test_0[1-5]*'    # a subset

# later, against a Python-port run with the same layout
python3 tests/reference/compare_outputs.py tests/reference/golden /path/to/python_run -q
# debugging only (not the standard): loosen floats and list order
python3 tests/reference/compare_outputs.py gold.tsv cand.tsv --unordered-lists --rel-tol 1e-9
```
