"""Compare the Python input reading with the original Tcl procs.

Runs only where ``tclsh`` is installed: the Tcl procs are sourced from
``share/tcl/AnnotSV`` and run on the test inputs of ``tests/AnnotSV``.
"""
import glob
import os
import shutil
import subprocess
import sys

import pytest

from annotsv.config import Config, ConfigExit
from annotsv.filtered_vcf import FilteredVCF
from annotsv.sv_input import add_na_in_samples_id_bed_col, create_bed_input_header_file
from annotsv.vcf import vcfs_to_bed
from annotsv import tcl

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
HERE = os.path.join(os.path.dirname(__file__), "tcl_parity")
TESTS = os.path.join(ROOT, "tests", "AnnotSV")
pytestmark = pytest.mark.skipif(shutil.which("tclsh") is None, reason="tclsh not installed")

VCFS = sorted(p for p in glob.glob(f"{TESTS}/**/*.vcf", recursive=True))
BEDS = sorted(p for p in glob.glob(f"{TESTS}/**/*.bed", recursive=True) if not p.endswith(".NA.bed"))


def _tclsh(script, *args):
    subprocess.run(["tclsh", os.path.join(HERE, script), ROOT, *map(str, args)],
                   check=True, capture_output=True, text=True)


def _collect(out_dir):
    """Sorted BED + unannotated lines (rescue lines come in hash order in Tcl)."""
    res = {}
    for name in sorted(os.listdir(out_dir)):
        path = os.path.join(out_dir, name)
        key = "bed" if name.endswith("_inputSVfile.bed") or name == "out.bed" else name
        if key in ("bed", "x.unannotated.tsv"):
            res[key] = sorted(open(path).read().split("\n"))
        else:
            res[key] = open(path).read()
    return res


@pytest.mark.parametrize("opts", [("1", "ALL", "1"), ("0", "DEL,INV", "0")])
@pytest.mark.parametrize("vcf", VCFS, ids=lambda p: os.path.relpath(p, TESTS))
def test_vcfs_to_bed_matches_tcl(tmp_path, vcf, opts):
    info, mode, ci = opts
    t, p = tmp_path / "t", tmp_path / "p"
    t.mkdir(), p.mkdir()
    _tclsh("vcfs_to_bed.tcl", vcf, t, info, 50, mode, ci)
    c = Config({"outputDir": str(p), "outputFile": "x.annotated.tsv", "SVminSize": "50",
                "bracketedAltMaxSize": "10", "bracketedAltMode": mode, "SVinputInfo": info,
                "includeCI": ci, "missingGTinSamplesid": "1", "candidateSnvIndelFiles": []})
    try:
        r = vcfs_to_bed(c, vcf)
        os.rename(r.bed_file, p / "out.bed")
        with open(p / "meta", "w") as f:
            f.write(r.vcf_header + "\n")
            for k in sorted(r.svlen):
                f.write(f"{k}={r.svlen[k]}\n")
    except ConfigExit as e:
        (p / "EXIT").write_text(f"{e.code}\n")
    assert _collect(p) == _collect(t)


@pytest.mark.parametrize("cols", [("-1", "-1"), ("4", "5")])
@pytest.mark.parametrize("bed", BEDS, ids=lambda p: os.path.relpath(p, TESTS))
def test_bed_header_matches_tcl(tmp_path, bed, cols):
    results = []
    for side in ("t", "p"):
        work, out = tmp_path / side / "in", tmp_path / side / "out"
        work.mkdir(parents=True), out.mkdir()
        shutil.copy(bed, work / "in.bed")
        if side == "t":
            _tclsh("bed_header.tcl", work / "in.bed", out, *cols)
        else:
            c = Config({"outputDir": str(out), "bedFile": str(work / "in.bed"),
                        "svtBEDcol": cols[0], "samplesidBEDcol": cols[1]})
            add_na_in_samples_id_bed_col(c, create_bed_input_header_file(c))
        results.append({n: open(os.path.join(d, n)).read()
                        for d in (work, out) for n in sorted(os.listdir(d))})
    assert results[0] == results[1]


@pytest.mark.parametrize("vcf,samples", [
    ("test_24_HtzCompound/input/3samples.snvindel.vcf", "prob father mother"),
    ("test_24_HtzCompound/input/test1.snvindel.vcf", "sample-001 sample-002"),
])
def test_filtered_vcf_matches_tcl(vcf, samples):
    vcf = os.path.join(TESTS, vcf)
    tcl_out = subprocess.run(["tclsh", os.path.join(HERE, "filtered_vcf.tcl"), ROOT, vcf, samples, "0"],
                             check=True, capture_output=True, text=True).stdout
    c = Config({"SVinputFile": "x.vcf", "candidateSnvIndelFiles": [vcf], "candidateSnvIndelSamples": samples.split(),
                "SVminSize": "50", "snvIndelPASS": "0", "formatColNumber": 1})
    fv = FilteredVCF(c)
    lines = [f"FULL={fv.annotate('FULL')}|"]
    key = lambda k: "FULL" if k == "FULL" else f"{k[0]},{k[1]}"  # noqa: E731
    lines += [f"{key(k)}={tcl.format_list(fv.htz_pos[k])}" for k in sorted(fv.htz_pos, key=key)]
    s = c["candidateSnvIndelSamples"]
    lines.append("samples=" + " ".join(s))
    hdr = "A\tFORMAT\t" + "\t".join(s) + "\tAnnotation_mode"
    line = ["x", "GT"] + ["0/1"] * len(s) + ["z"] * 13
    for ch, a, b in [("1", 1, 300000000), ("1", 10000, 20000000), ("2", 1, 100000000), ("X", 1, 200000000)]:
        lines.append(f"{ch}:{a}-{b}={fv.annotate(ch, a, b, line, hdr)}|")
    expected = [l for l in tcl_out.split("\n") if l and not l.startswith(("...", "\t"))]
    assert "\n".join(lines).split("\n") == expected
