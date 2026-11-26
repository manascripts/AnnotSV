#!/usr/bin/env bash
# Run the Tcl AnnotSV test suite (tests/AnnotSV/test_*/command_Tcl.sh) on a
# scratch copy of the test inputs and collect every output file as the golden
# reference for the Python port.
#
# Usage: tests/reference/run_tcl_reference.sh [test_dir_glob ...]
#   e.g. tests/reference/run_tcl_reference.sh            # all tests
#        tests/reference/run_tcl_reference.sh 'test_0*'  # a subset
#
# Env:
#   ANNOTSV   AnnotSV install dir (default: repo root)
#   WORKDIR   scratch dir for the copied tests (default: mktemp)
#   GOLDEN    destination of the golden outputs (default: tests/reference/golden)

set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export ANNOTSV="${ANNOTSV:-$REPO}"
GOLDEN="${GOLDEN:-$REPO/tests/reference/golden}"
WORKDIR="${WORKDIR:-$(mktemp -d -t annotsv_ref.XXXXXX)}"

patterns=("$@")
[ ${#patterns[@]} -eq 0 ] && patterns=("test_*")

# The test scripts modify their own inputs (gzip, temp configfiles), so work on a copy.
cp -r "$REPO/tests/AnnotSV/." "$WORKDIR/"
mkdir -p "$GOLDEN"
summary="$GOLDEN/SUMMARY.tsv"
[ -f "$summary" ] || printf "test\tstatus\texit_code\tseconds\n" > "$summary"

for pat in "${patterns[@]}"; do
  for d in "$WORKDIR"/$pat/; do
    [ -f "$d/command_Tcl.sh" ] || continue
    t="$(basename "$d")"
    echo "######## $t  ($(date +%T))"
    touch "$d/.run_marker"; sleep 1
    start=$(date +%s)
    ( cd "$d" && bash ./command_Tcl.sh ) > "$d/command_Tcl.log" 2>&1
    rc=$?
    secs=$(( $(date +%s) - start ))
    if [ $rc -eq 0 ] && grep -q "ok - Finished" "$d/command_Tcl.log"; then status=ok; else status=FAIL; fi
    echo "   -> $status (exit $rc, ${secs}s)"

    dest="$GOLDEN/$t"
    rm -rf "$dest"; mkdir -p "$dest"
    # Every file the run produced: ./output/*, plus anything written beside the inputs
    # (some tests write *.annotated.tsv next to command_Tcl.sh). Timestamped AnnotSV temp files are skipped.
    ( cd "$d" && find . -type f -newer .run_marker ! -name command_Tcl.log ! -name .run_marker ! -name '[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]-[0-9]*_AnnotSV_*' -print0 ) |
      while IFS= read -r -d '' f; do
        mkdir -p "$dest/$(dirname "$f")"; cp -p "$d/$f" "$dest/$f"
      done
    # Log with absolute scratch paths normalised so it diffs cleanly between runs.
    sed -e "s#$WORKDIR#<WORKDIR>#g" -e "s#$ANNOTSV#<ANNOTSV>#g" "$d/command_Tcl.log" > "$dest/command_Tcl.log"

    grep -v "^$t	" "$summary" > "$summary.tmp" || true
    printf "%s\t%s\t%s\t%s\n" "$t" "$status" "$rc" "$secs" >> "$summary.tmp"
    { head -1 "$summary.tmp"; tail -n +2 "$summary.tmp" | sort; } > "$summary"; rm -f "$summary.tmp"
  done
done

echo "Golden outputs in $GOLDEN (scratch copy: $WORKDIR)"
