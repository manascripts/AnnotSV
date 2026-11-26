"""Compare the Python region annotations with the Tcl procs.

Covers checkCytoband, checkRepeatFile, checkSegDupFile, checkGapFile,
checkENCODEblacklistFile, checkFASTAfiles (GC content), checkTADfiles,
checkCOSMICfile, checkUsersBED and the matching columns of OrganizeAnnotation
(CytoBand, user SVincludedInFt / FtIncludedInSV / AnyOverlap / BNDproximity
BED files, TAD, GC content, Repeat, Gap, SegDup, ENCODE blacklist).

Runs only where ``tclsh`` and ``bedtools`` are installed. Every download is
synthetic, built around the breakpoints of each input.
"""
import gzip
import io
import os
import random
import re
import shutil
import subprocess
import tarfile

import pytest

from test_genes_parity import INPUTS, TESTS, ROOT, _annotations as _genes_annotations, _breakpoints, _collect
from test_gene_based_parity import _tree
from annotsv import regions
from annotsv.closest_genes import ClosestGenes
from annotsv.config import Config, ConfigExit
from annotsv.genes import check_genes_refseq_file, genes_annotation
from annotsv.sv_input import read_sv_input
from annotsv.write import organize_annotation

DRIVER = os.path.join(os.path.dirname(__file__), "tcl_parity", "regions.tcl")
pytestmark = pytest.mark.skipif(shutil.which("tclsh") is None or shutil.which("bedtools") is None,
                                reason="tclsh or bedtools not installed")
SAMPLE = [p for p in INPUTS if re.search(r"test_(02_HG00096/input/.*\.bed|16_reciprocalYes/input/test\.vcf|"
                                         r"22_knotAnnotSV_41SV/input/test\.41)", p)]
ALL_COLS = ["CytoBand", "TAD_coordinate", "ENCODE_experiment", "GC_content_left", "GC_content_right",
            "Repeat_coord_left", "Repeat_type_left", "Repeat_coord_right", "Repeat_type_right", "Gap_left",
            "Gap_right", "SegDup_left", "SegDup_right", "ENCODE_blacklist_left",
            "ENCODE_blacklist_characteristics_left", "ENCODE_blacklist_right",
            "ENCODE_blacklist_characteristics_right"]
SOME_COLS = ["ENCODE_experiment", "GC_content_right", "Gap_left", "ENCODE_blacklist_right"]
FA_MAX = 300_000


def _small_input(path):
    """An SV BED input with small coordinates, so that the GC content is computed."""
    rnd = random.Random("small")
    lines = ["#chrom\tstart\tend\tid\tSV_type"]
    for k in range(40):
        chrom = rnd.choice(["chr1", "chr2", "chrX", "3"])
        start = rnd.randint(0, FA_MAX - 2000)
        end = start + rnd.choice([1, 60, 400, 5000, 80000])
        lines.append(f"{chrom}\t{start}\t{min(end, FA_MAX + 300)}\tsv{k}\t{rnd.choice(['DEL', 'DUP', 'INS', 'INV'])}")
    lines.append(f"chr2\t{FA_MAX - 120}\t{FA_MAX + 50}\tsvEnd\tDEL")
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def _gz(path):
    return io.TextIOWrapper(gzip.GzipFile(path, "wb", mtime=0))


def _region(rnd, pos, small=50, big=60000):
    start = max(0, pos + rnd.randint(-big, big // 4))
    return start, start + rnd.randint(small, big)


def _annotations(root, input_file):
    rnd = random.Random("regions" + input_file)
    ann, _ = _genes_annotations(root, input_file, "RefSeq")
    base = root / "ann" / "Annotations_Human"
    bps = _breakpoints(input_file)
    chroms = sorted({re.sub("^chr", "", c) for c, _ in bps})

    def bp_file(name, n_extra, row):
        d = base / "BreakpointsAnnotations" / name / "GRCh38"
        d.mkdir(parents=True)
        lines = ["#chrom\tstart\tend\tname"]
        for chrom, pos in bps:
            for _ in range(rnd.choice([0, 0, 1, 1, 2])):
                s, e = _region(rnd, pos, 10, 300)
                lines.append(row("chr" + re.sub("^chr", "", chrom), s, e))
        lines += [row(rnd.choice(["chrUn_x", "chr23", "chr1_alt"]), 10, 500) for _ in range(n_extra)]
        rnd.shuffle(lines[1:])
        (d / f"{name}.bed").write_text("\n".join(lines) + "\n\n")

    bp_file("Gap", 2, lambda c, s, e: f"{c}\t{s}\t{e}\t{rnd.choice(['contig', 'telomere'])}")
    bp_file("SegDup", 1, lambda c, s, e: f"{c}\t{s}\t{e}\tsd{rnd.randint(0, 99)}")
    bp_file("Repeat", 3, lambda c, s, e: f"{c}\t{s}\t{e}\t{rnd.choice(['L1', 'Alu', '(CA)n', 'MIR'])}")
    bp_file("ENCODEblacklist", 1, lambda c, s, e: f"{c}\t{s}\t{e}\t{rnd.choice(['High Signal Region', 'Low Mappability'])}")

    # Cytoband
    d = base / "AnyOverlap" / "CytoBand" / "GRCh38"
    d.mkdir(parents=True)
    cyto = []
    for chrom in chroms:
        top = max(p for c, p in bps if re.sub("^chr", "", c) == chrom) + 5_000_000
        s, k = 0, 0
        while s < top:
            e = s + rnd.randint(200_000, 3_000_000)
            cyto.append(f"chr{chrom}\t{s}\t{e}\t{'p' if k < 5 else 'q'}{k}.{rnd.randint(1, 3)}")
            s, k = e, k + 1
    (d / "cytoBand_GRCh38.bed").write_text("\n".join(cyto) + "\n")

    # TAD
    d = base / "FtIncludedInSV" / "TAD" / "GRCh38"
    d.mkdir(parents=True)
    for exp in ("ENCFF001AAA", "ENCFF002BBB"):
        with _gz(d / f"{exp}.bed.gz") as f:
            f.write("#tad\n")
            for chrom, pos in bps:
                for _ in range(rnd.choice([0, 1, 2])):
                    s = max(0, pos + rnd.randint(-100000, 100000))
                    name = rnd.choice(["tad", "x___boundary"])
                    f.write(f"chr{re.sub('^chr', '', chrom)}\t{s}\t{s + rnd.randint(1000, 400000)}\t{name}\n")
            f.write("chr1\t10\t20000\ttad\n\n")

    # COSMIC
    d = base / "FtIncludedInSV" / "COSMIC" / "GRCh38"
    d.mkdir(parents=True)
    (base / "Users" / "GRCh38" / "FtIncludedInSV").mkdir(parents=True)
    old_format = input_file.endswith("small.bed")   # "Chromosome:G_Start..G_Stop" (before 2023)
    with _gz(d / "CosmicCompleteCNA.tsv.gz") as f:
        if old_format:
            f.write("CNV_ID\tSAMPLE\tMUT_TYPE\tChromosome:G_Start..G_Stop\n")
        else:
            f.write("COSMIC_CNV_ID\tSAMPLE\tMUT_TYPE\tCHROMOSOME\tGENOME_START\tGENOME_STOP\n")
        for k, (chrom, pos) in enumerate(bps[:80]):
            s, e = _region(rnd, pos, 100, 3000)
            chrom = re.sub('^chr', '', chrom)
            if old_format:
                coord = f"{chrom}:{s}..{e}" if k % 9 != 1 else "bad"   # Tcl keeps the previous coordinates
                row = f"COSV{k % 13}\tS\t{rnd.choice(['gain', 'loss'])}\t{coord}\n"
            else:
                row = f"COSV{k % 13}\tS\t{rnd.choice(['gain', 'loss'])}\t{chrom}\t{s}\t{e}\n"
            f.write(row + (row if rnd.random() < 0.2 else ""))

    # Users BED files
    users = base / "Users" / "GRCh38"

    def user_file(kind, name, header, ncols, big):
        d = users / kind
        d.mkdir(exist_ok=True)
        lines = [header] if header else []
        for chrom, pos in bps:
            for _ in range(rnd.choice([0, 1, 1, 2])):
                s, e = _region(rnd, pos, 1, big)
                extra = [rnd.choice(["a", "b b", "z", f"v{rnd.randint(0, 5)}"]) for _ in range(ncols - 3)]
                lines.append("\t".join([f"chr{re.sub('^chr', '', chrom)}", str(s), str(e)] + extra))
        # Tcl crashes on a user BED without data lines (lrepeat on an empty list in checkBed)
        lines.append("\t".join(["chr22", "100", "200"] + ["w"] * (ncols - 3)))
        (d / f"{name}.bed").write_text("\n".join(lines) + "\n")

    user_file("SVincludedInFt", "svin", "#chr\tstart\tend\tsvinA\tsvinB", 5, 2_000_000)
    user_file("FtIncludedInSV", "ftin", "", 4, 5000)
    user_file("AnyOverlap", "any", "#c\ts\te\tanyA", 4, 30000)
    user_file("BNDproximity", "bnd", "#c\ts\te\tbndA\tbndB", 5, 150)
    user_file("AnyOverlap", "bad", "#c\ts\te\tbadA", 4, 30000)
    (users / "AnyOverlap" / "bad.header.tsv").write_text("x\ty\n")

    # GC content
    d = base / "BreakpointsAnnotations" / "GCcontent" / "GRCh38"
    d.mkdir(parents=True)
    with tarfile.open(d / "hg38.chromFa.tar.gz", "w:gz", format=tarfile.USTAR_FORMAT) as tar:
        for chrom in chroms + (["Y"] if "Y" not in chroms else []):
            seq = "".join(rnd.choice("ACGTacgtNNR") for _ in range(FA_MAX))
            text = f">chr{chrom}\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n"
            data = text.encode()
            info = tarfile.TarInfo(f"chroms/chr{chrom}.fa")
            info.size, info.mtime = len(data), 0
            tar.addfile(info, io.BytesIO(data))
    return ann


def _python(ann, out, input_file, mode, metrics, overlap, reciprocal, cols):
    c = Config({"annotationsDir": ann, "organism": "Human", "genomeBuild": "GRCh38", "outputDir": out,
                "outputFile": "x.annotated.tsv", "SVinputFile": input_file, "svtBEDcol": "4",
                "samplesidBEDcol": "-1", "txFile": "", "annotationMode": mode, "tx": "RefSeq",
                "SVinputInfo": "1", "SVminSize": "50", "bracketedAltMaxSize": "10", "bracketedAltMode": "ALL",
                "includeCI": "1", "missingGTinSamplesid": "1", "candidateSnvIndelFiles": [], "metrics": metrics,
                "overlap": overlap, "reciprocal": reciprocal, "breakpointProximity": "100", "extann": [],
                "bedtools": "bedtools", "outputColHeader": ["AnnotSV_ID", "Samples_ID"] + cols})
    try:
        sv = read_sv_input(c)
        check_genes_refseq_file(c)
        regions.check_fasta_files(c)
        regions.check_breakpoint_files(c)
        regions.check_cytoband(c)
        regions.check_tad_files(c)
        regions.check_cosmic_file(c)
        regions.check_users_bed(c)
        g = genes_annotation(c)
        organize_annotation(c, g, vcf_header=sv.vcf_header, svlen=sv.svlen, closest_genes=ClosestGenes(c),
                            header_file_to_remove=sv.header_file_to_remove,
                            regions=regions.RegionAnnotations(c))
    except ConfigExit as e:
        with open(os.path.join(out, "EXIT"), "w") as f:
            f.write(f"{e.code}\n")


CASES = [  # annotationMode metrics overlap reciprocal columns
    ("both", "us", "100", "0", "all"),
    ("split", "fr", "70", "1", "all"),
    ("full", "us", "50", "0", "some"),
]


@pytest.fixture(scope="module")
def small_input(tmp_path_factory):
    return _small_input(tmp_path_factory.mktemp("small") / "small.bed")


@pytest.mark.parametrize("case", CASES, ids="-".join)
@pytest.mark.parametrize("which", ["small"] + SAMPLE, ids=lambda p: p if p == "small" else os.path.relpath(p, TESTS))
def test_regions_match_tcl(tmp_path, which, case, small_input):
    input_file = small_input if which == "small" else which
    mode, metrics, overlap, reciprocal, colset = case
    cols = ALL_COLS if colset == "all" else SOME_COLS
    results = []
    for side in ("t", "p"):
        root = tmp_path / side
        root.mkdir()
        ann = _annotations(root, input_file)
        work, out = root / "in", root / "out"
        work.mkdir(), out.mkdir()
        inp = str(work / os.path.basename(input_file))
        shutil.copy(input_file, inp)
        if side == "t":
            r = subprocess.run(["tclsh", DRIVER, ROOT, ann, str(out), inp, "GRCh38", mode, metrics, overlap,
                                reciprocal, " ".join(cols)], capture_output=True, text=True)
            assert r.returncode == 0 and "Error" not in r.stdout[-300:], (r.stdout + r.stderr)[-3000:]
            tsv = out / "x.annotated.tsv"
            if tsv.exists():
                tsv.write_text(re.sub("\t@STUB@", "", tsv.read_text()))
        else:
            _python(ann, str(out), inp, mode, metrics, overlap, reciprocal, cols)
        results.append(root)
    out_p = results[1] / "out" / "x.annotated.tsv"
    header = open(out_p).readline() if out_p.exists() else ""
    # The ranking columns come after the columns not ported in this test: checked in test_ranking_parity
    n_cols = len(header.rstrip("\n").split("\t")) - (3 if header.rstrip("\n").endswith("ACMG_class") else 0)
    t, p = (_collect(r / "out", n_cols) for r in results)
    assert sorted(t) == sorted(p)
    for name in t:
        assert t[name].split("\n") == p[name].split("\n"), name
    t_files, p_files = _tree(results[0]), _tree(results[1])
    assert sorted(t_files) == sorted(p_files)
    for rel in t_files:
        assert t_files[rel] == p_files[rel], rel
    assert "Gap_left" in header and "bndA" in header
