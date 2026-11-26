"""Compare the Python gene annotation and core output columns with the Tcl
procs (checkGenes*File, genesAnnotation, OrganizeAnnotation).

Runs only where ``tclsh`` and ``bedtools`` are installed. The annotation
bundle is not needed: a synthetic UCSC RefSeq table is generated around the
breakpoints of each test input (several transcripts per gene, both strands,
non-coding transcripts, duplicates, non-standard contigs, versions).
"""
import glob
import gzip
import os
import random
import re
import shutil
import subprocess

import pytest

from annotsv.closest_genes import ClosestGenes
from annotsv.config import Config, ConfigExit
from annotsv.genes import check_genes_file, genes_annotation
from annotsv.sv_input import read_sv_input
from annotsv.write import organize_annotation

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DRIVER = os.path.join(os.path.dirname(__file__), "tcl_parity", "genes_annotation.tcl")
TESTS = os.path.join(ROOT, "tests", "AnnotSV")
pytestmark = pytest.mark.skipif(shutil.which("tclsh") is None or shutil.which("bedtools") is None,
                                reason="tclsh or bedtools not installed")

INPUTS = sorted(p for p in glob.glob(f"{TESTS}/**/*.vcf", recursive=True) + glob.glob(f"{TESTS}/**/*.bed", recursive=True)
                if not p.endswith(".NA.bed") and os.path.getsize(p) < 2_000_000)


def _breakpoints(path):
    out = []
    with open(path, errors="replace") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            ls = re.split(r"\t| +", line.rstrip("\n")) if path.endswith(".bed") else line.split("\t")
            if len(ls) < 3:
                continue
            chrom = ls[0]
            try:
                pos = [int(ls[1])]
                if path.endswith(".bed"):
                    pos.append(int(ls[2]))
                else:
                    m = re.search(r"(?:^|;)END=(\d+)", ls[7] if len(ls) > 7 else "")
                    if m:
                        pos.append(int(m.group(1)))
            except ValueError:
                continue
            out += [(chrom, p) for p in pos]
    return out[:400]


def _transcript(rnd, name, chrom, strand, gene, center):
    n_ex = rnd.randint(1, 6)
    start = max(0, center - rnd.randint(0, 40000))
    bounds, pos = [], start
    for _ in range(n_ex):
        ex = rnd.randint(50, 3000)
        bounds.append((pos, pos + ex))
        pos += ex + rnd.randint(100, 15000)
    tx_start, tx_end = bounds[0][0], bounds[-1][1]
    if rnd.random() < 0.25:
        cds_start = cds_end = tx_end                       # non-coding (NR_)
    else:
        cds_start = rnd.randint(tx_start, bounds[0][1])
        cds_end = rnd.randint(bounds[-1][0], tx_end)
    starts = "".join(f"{a}," for a, _ in bounds)
    ends = "".join(f"{b}," for _, b in bounds)
    return [str(rnd.randint(0, 2000)), name, chrom, strand, str(tx_start), str(tx_end), str(cds_start), str(cds_end),
            str(n_ex), starts, ends, "0", gene, "cmpl", "cmpl", ""]


def _genes_rows(input_file):
    rnd = random.Random(input_file)
    rows, k = [], 0
    for chrom, pos in _breakpoints(input_file):
        chrom = "chr" + re.sub("^chr", "", chrom)
        for _ in range(rnd.randint(1, 2)):
            gene, strand = f"G{k}", rnd.choice("+-")
            k += 1
            for t in range(rnd.randint(1, 3)):
                name = f"NM_{k:05d}{t}.{rnd.randint(1, 9)}" if rnd.random() < 0.9 else f"NR_{k:05d}{t}"
                rows.append(_transcript(rnd, name, chrom, strand, gene, pos + rnd.randint(-5000, 5000)))
        for _ in range(rnd.randint(0, 2)):             # nearest genes, inside or beyond 5 Mb
            gene = f"F{k}" if rnd.random() < 0.8 else f"G{max(0, k - 1)}"
            k += 1
            rows.append(_transcript(rnd, f"NM_8{k:05d}.1", chrom, rnd.choice("+-"), gene,
                                    max(0, pos + rnd.choice((-1, 1)) * rnd.randint(20000, 6000000))))
        if rnd.random() < 0.1:
            rows.append(list(rows[-1]))                    # duplicated line
            rows.append(_transcript(rnd, f"NM_9{k}.1", chrom + "_alt", "+", "ALT", pos))
    return rows


def _annotations(tmp_path, input_file, tx):
    gdir = tmp_path / "ann" / "Annotations_Human" / "Genes" / "GRCh38"
    gdir.mkdir(parents=True)
    rows = _genes_rows(input_file)
    if tx == "RefSeq":
        with gzip.open(gdir / "ncbiRefSeq.txt.gz", "wt") as f:
            f.write("".join("\t".join(r) + "\n" for r in rows))
    else:
        bed = [[r[2][3:], r[4], r[5], r[3], r[12], r[1].replace("NM_", "ENST"), r[6], r[7], r[9], r[10]]
               for r in rows if not r[2].endswith("_alt")]
        bed.sort(key=lambda r: (r[0], int(r[1])))
        (gdir / "genes.ENSEMBL.sorted.tmp.bed").write_text("".join("\t".join(r) + "\n" for r in bed))
    return str(tmp_path / "ann"), [r[1] for r in rows]


def _python(ann, out, input_file, svt, sid, tx_file, mode, tx, info):
    c = Config({"annotationsDir": ann, "organism": "Human", "genomeBuild": "GRCh38", "outputDir": out,
                "outputFile": "x.annotated.tsv", "SVinputFile": input_file, "svtBEDcol": svt,
                "samplesidBEDcol": sid, "txFile": tx_file, "annotationMode": mode, "tx": tx,
                "SVinputInfo": info, "SVminSize": "50", "bracketedAltMaxSize": "10", "bracketedAltMode": "ALL",
                "includeCI": "1", "missingGTinSamplesid": "1", "candidateSnvIndelFiles": [],
                "outputColHeader": ["AnnotSV_ID", "Samples_ID"]})
    try:
        sv = read_sv_input(c)
        check_genes_file(c)
        g = genes_annotation(c)
        shutil.copy(g.full_and_split_file, os.path.join(out, "fullsplit.txt"))
        with open(os.path.join(out, "genes.txt"), "w") as f:
            for k in sorted(g.genes_by_sv):
                f.write(f"{k}={g.genes_by_sv[k]}\n")
        organize_annotation(c, g, vcf_header=sv.vcf_header, svlen=sv.svlen, closest_genes=ClosestGenes(c),
                            header_file_to_remove=sv.header_file_to_remove)
    except ConfigExit as e:
        with open(os.path.join(out, "EXIT"), "w") as f:
            f.write(f"{e.code}\n")


def _collect(out_dir, n_cols):
    res = {}
    for name in sorted(os.listdir(out_dir)):
        path = os.path.join(out_dir, name)
        if name.endswith(".bash") or os.path.isdir(path):
            continue
        text = open(path).read()
        if name == "x.annotated.tsv":       # Tcl also writes the not yet ported columns
            text = "\n".join("\t".join(l.split("\t")[:n_cols]) for l in text.split("\n"))
        if name.endswith(".unannotated.tsv"):    # reciprocal breakends come in Tcl hash order
            text = "\n".join(sorted(text.split("\n")))
        res[re.sub(r"^[0-9]{8}-[0-9]{6}_[0-9]+_", "", name)] = text   # timestamp and pid
    return res


CASES = [  # svtBEDcol samplesidBEDcol txFile annotationMode tx SVinputInfo
    ("4", "-1", False, "both", "RefSeq", "1"),
    ("-1", "-1", True, "split", "RefSeq", "0"),
    ("4", "-1", False, "full", "ENSEMBL", "0"),
]


@pytest.mark.parametrize("case", CASES, ids=lambda c: "-".join(x if isinstance(x, str) else str(x) for x in c))
@pytest.mark.parametrize("input_file", INPUTS, ids=lambda p: os.path.relpath(p, TESTS))
def test_genes_annotation_matches_tcl(tmp_path, input_file, case):
    svt, sid, use_tx_file, mode, tx, info = case
    results = []
    for side in ("t", "p"):
        root = tmp_path / side
        root.mkdir()
        ann, names = _annotations(root, input_file, tx)
        tx_file = ""
        if use_tx_file:
            tx_file = str(root / "tx.txt")
            with open(tx_file, "w") as f:
                f.write(" ".join(n for n in names[::3]) + "\n" + "|".join(names[1::7]) + "\n")
        work, out = root / "in", root / "out"
        work.mkdir(), out.mkdir()
        inp = str(work / os.path.basename(input_file))
        shutil.copy(input_file, inp)
        if side == "t":
            r = subprocess.run(["tclsh", DRIVER, ROOT, ann, str(out), inp, svt, sid, tx_file, mode, tx, info],
                               capture_output=True, text=True)
            assert r.returncode == 0 and "Error" not in r.stdout[-300:], (r.stdout + r.stderr)[-3000:]
        else:
            _python(ann, str(out), inp, svt, sid, tx_file, mode, tx, info)
        results.append(root)
    header = open(results[1] / "out" / "x.annotated.tsv").readline() if (results[1] / "out" / "x.annotated.tsv").exists() else ""
    # The ranking columns come after the columns not ported in this test: checked in test_ranking_parity
    n_cols = len(header.rstrip("\n").split("\t")) - (3 if header.rstrip("\n").endswith("ACMG_class") else 0)
    t, p = (_collect(r / "out", n_cols) for r in results)
    assert sorted(t) == sorted(p)
    for name in t:
        assert t[name].split("\n") == p[name].split("\n"), name
    for sub in ("in", "ann"):
        for d, _, files in os.walk(results[0] / sub):
            for fn in files:
                rel = os.path.relpath(os.path.join(d, fn), results[0])
                a, b = open(os.path.join(d, fn), "rb").read(), open(results[1] / rel, "rb").read()
                assert a == b, rel
