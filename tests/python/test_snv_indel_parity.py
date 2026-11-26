"""Compare the Python SNV/indel counts and compound htz columns with the Tcl procs.

Covers VCFannotation (Count_hom / Count_htz / Count_htz/allHom /
Count_htz/total(cohort) per sample, Count_total(cohort); bcftools norm of the
multiallelic sites, "chr" removal, bedtools intersect with the SV BED, the SV,
FILTER and GT filters, the >50 variants and deletion-only rules) and the
compound_htz columns written by OrganizeAnnotation through
filteredVCFannotation (full and split lines).

Runs only where ``tclsh``, ``bedtools`` and ``bcftools`` are installed. The
SNV/indel VCFs are synthetic, built around the breakpoints of each input.
"""
import gzip
import os
import random
import re
import shutil
import subprocess

import pytest

from test_genes_parity import TESTS, ROOT, _annotations as _genes_annotations, _collect
from test_gene_based_parity import _tree
from annotsv.closest_genes import ClosestGenes
from annotsv.config import Config, ConfigExit
from annotsv.filtered_vcf import FilteredVCF
from annotsv.genes import check_genes_refseq_file, genes_annotation
from annotsv.sv_input import read_sv_input
from annotsv.write import organize_annotation

DRIVER = os.path.join(os.path.dirname(__file__), "tcl_parity", "snv_indel.tcl")
pytestmark = pytest.mark.skipif(any(shutil.which(t) is None for t in ("tclsh", "bedtools", "bcftools")),
                                reason="tclsh, bedtools or bcftools not installed")

INPUTS = [  # input, svtBEDcol
    ("test_24_HtzCompound/input/3samples.severalSV.vcf", "4"),
    ("test_24_HtzCompound/input/test1.SV.vcf", "4"),
    ("test_26_HtzHomCounts_deletionFiltering/input/3samples.SV.vcf", "4"),
    ("test_26_HtzHomCounts_deletionFiltering/input/SV.bed", "3"),
]
CASES = [("both", "0"), ("split", "1"), ("full", "0")]  # annotationMode snvIndelPASS
GTS = ["0/1", "1/1", "0|1", "1|0", "1|1", "0/0", "./.", "0/.", "1", ".", "0"]


def _svs(path):
    """(chrom, start, end) of each SV of the input."""
    out = []
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        ls = line.rstrip("\n").split("\t")
        if path.endswith(".bed"):
            out.append((ls[0], int(ls[1]), int(ls[2])))
        else:
            m = re.search(r"(?:^|;)END=(\d+)", ls[7])
            out.append((ls[0], int(ls[1]), int(m.group(1)) if m else int(ls[1]) + 1))
    return out


def _samples(path):
    if path.endswith(".bed"):
        return ["s1", "s2"]
    for line in open(path):
        if line.startswith("#CHROM"):
            return line.rstrip("\n").split("\t")[9:]
    return []


def _positions(rnd, svs):
    pos = set()
    for chrom, start, end in svs:
        for p in (start - 1, start, start + 1, end - 1, end, end + 1):
            pos.add((chrom, p))
        for _ in range(rnd.randint(0, 8)):
            pos.add((chrom, max(1, start + rnd.randint(-2000, 2000))))
        if end - start < 3_000_000:
            for _ in range(rnd.choice([3, 20, 70, 120])):
                pos.add((chrom, rnd.randint(start, end)))
    return pos


def _row(rnd, chrom, pos, n_samples):
    kind = rnd.random()
    ref, alt, fmt, info = "A", rnd.choice(["C", "G", "T"]), "GT:DP", "."
    if kind < 0.08:
        alt = "C,G"                                         # multiallelic (bcftools norm)
    elif kind < 0.12:
        alt, info = "<DEL>", f"SVTYPE=DEL;END={pos + 300}"
    elif kind < 0.14:
        alt, info = "A]17:1584563]", "SVTYPE=BND"
    elif kind < 0.17:
        ref = "A" + "C" * rnd.choice([49, 50, 51, 70])      # indel / SV size
    elif kind < 0.2:
        alt = "AT"
    elif kind < 0.22:
        fmt = "DP"                                           # no GT
    filt = "PASS" if rnd.random() < 0.8 else rnd.choice([".", "LowQual"])
    gts = []
    for _ in range(n_samples):
        gt = rnd.choice(GTS)
        if "," in alt and rnd.random() < 0.5:
            gt = rnd.choice(["1/2", "0/2", "2/2"])
        gts.append(f"{gt}:7" if fmt == "GT:DP" else "7")
    return [chrom, str(pos), ".", ref, alt, "50", filt, info, fmt] + gts


def _vcf(path, rows, samples, gz):
    text = ("##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n"
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"DP\">\n"
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"END\">\n"
            "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SVTYPE\">\n"
            + "".join(f"##contig=<ID={c}>\n" for c in sorted({r[0] for r in rows} | {"17"})) +
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n"
            + "".join("\t".join(r) + "\n" for r in rows))
    if gz:
        with gzip.GzipFile(path, "wb", mtime=0) as f:
            f.write(text.encode())
    else:
        with open(path, "w") as f:
            f.write(text)


def _snv_files(work, input_file):
    rnd = random.Random("snvindel" + input_file)
    svs = _svs(input_file)
    samples = _samples(input_file)
    files = []
    # File 1: every sample of the SV input, sorted positions
    pos = sorted(_positions(rnd, svs), key=lambda x: (x[0], x[1]))
    rows = [_row(rnd, c, p, len(samples)) for c, p in pos]
    files.append((work / "a.snv.vcf", rows, samples, False))
    # File 2: one sample of the input and another one, "chr" prefix toggled, a few lines out of order
    s2 = samples[-1:] + ["other"]
    pos = sorted(_positions(rnd, svs), key=lambda x: (x[0], x[1]))
    rows = [_row(rnd, ("chr" + c) if not c.startswith("chr") else c[3:], p, len(s2)) for c, p in pos]
    if len(rows) > 10:
        rows[3], rows[7] = rows[7], rows[3]
    files.append((work / "b.snv.vcf.gz", rows, s2, True))
    for path, rows, s, gz in files:
        _vcf(path, rows, s, gz)
    all_samples = sorted(set(samples) | set(s2))
    return [str(f[0]) for f in files], all_samples


def _python(ann, out, inp, mode, svt, snv, samples, cand, pass_only):
    c = Config({"annotationsDir": ann, "organism": "Human", "genomeBuild": "GRCh38", "outputDir": out,
                "outputFile": "x.annotated.tsv", "SVinputFile": inp, "svtBEDcol": svt,
                "samplesidBEDcol": "-1", "txFile": "", "annotationMode": mode, "tx": "RefSeq",
                "SVinputInfo": "1", "SVminSize": "50", "bracketedAltMaxSize": "10", "bracketedAltMode": "ALL",
                "includeCI": "1", "missingGTinSamplesid": "1", "candidateSnvIndelFiles": cand,
                "candidateSnvIndelSamples": samples if cand else [], "snvIndelFiles": snv,
                "snvIndelSamples": samples, "snvIndelPASS": pass_only, "bedtools": "bedtools",
                "bcftools": "bcftools", "extann": [], "outputColHeader": ["AnnotSV_ID", "Samples_ID"]})
    try:
        if cand:
            c.filtered_vcf = FilteredVCF(c)
            c.filtered_vcf.annotate("FULL")
        sv = read_sv_input(c)
        check_genes_refseq_file(c)
        g = genes_annotation(c)
        organize_annotation(c, g, vcf_header=sv.vcf_header, svlen=sv.svlen, closest_genes=ClosestGenes(c),
                            header_file_to_remove=sv.header_file_to_remove)
    except ConfigExit as e:
        with open(os.path.join(out, "EXIT"), "w") as f:
            f.write(f"{e.code}\n")


@pytest.mark.parametrize("case", CASES, ids="-".join)
@pytest.mark.parametrize("which", INPUTS, ids=lambda w: w[0])
def test_snv_indel_match_tcl(tmp_path, which, case):
    rel, svt = which
    input_file = os.path.join(TESTS, rel)
    mode, pass_only = case
    results = []
    for side in ("t", "p"):
        root = tmp_path / side
        root.mkdir()
        ann, _ = _genes_annotations(root, input_file, "RefSeq")
        work, out = root / "in", root / "out"
        work.mkdir(), out.mkdir()
        inp = str(work / os.path.basename(input_file))
        shutil.copy(input_file, inp)
        snv, samples = _snv_files(work, input_file)
        cand = [] if inp.endswith(".bed") else snv
        if side == "t":
            r = subprocess.run(["tclsh", DRIVER, ROOT, ann, str(out), inp, mode, svt, " ".join(snv),
                                " ".join(samples), " ".join(cand), " ".join(samples) if cand else "",
                                pass_only], capture_output=True, text=True)
            assert r.returncode == 0 and "Error" not in r.stdout[-300:], (r.stdout + r.stderr)[-3000:]
            tsv = out / "x.annotated.tsv"
            if tsv.exists():
                tsv.write_text(re.sub("\t@STUB@", "", tsv.read_text()))
        else:
            _python(ann, str(out), inp, mode, svt, snv, samples, cand, pass_only)
        results.append(root)
    out_p = results[1] / "out" / "x.annotated.tsv"
    header = open(out_p).readline() if out_p.exists() else ""
    n_cols = len(header.rstrip("\n").split("\t"))
    t, p = (_collect(r / "out", n_cols) for r in results)
    assert sorted(t) == sorted(p)
    for name in t:
        assert t[name].split("\n") == p[name].split("\n"), name
    t_files, p_files = _tree(results[0]), _tree(results[1])
    assert sorted(t_files) == sorted(p_files)
    for rel_path in t_files:
        assert t_files[rel_path] == p_files[rel_path], rel_path
    assert "Count_total(cohort)" in header
    if cand:
        assert "compound_htz(" in header
