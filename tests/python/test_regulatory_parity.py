"""Compare the Python regulatory elements annotation (check*files,
regulatoryElementsAnnotation, RE_gene) with the Tcl procs.

Runs only where ``tclsh`` and ``bedtools`` are installed. The downloads are
synthetic: EnhancerAtlas, GeneHancer, miRTargetLink, ABC (raw and an unsorted
formatted file), MPRA, ClinGen, OMIM morbid genes and a candidate genes file,
built around the breakpoints of each test input with the genes of the
synthetic genes table (plus unknown genes and badly formatted lines).
"""
import gzip
import io
import os
import random
import re
import shutil
import subprocess

import pytest

from test_genes_parity import INPUTS, TESTS, ROOT, _breakpoints, _collect, _genes_rows
from annotsv.closest_genes import ClosestGenes
from annotsv.config import Config, ConfigExit
from annotsv.genes import check_genes_ensembl_file, check_genes_refseq_file, genes_annotation
from annotsv import regulatory_elements as rem
from annotsv.sv_input import read_sv_input
from annotsv.write import organize_annotation

DRIVER = os.path.join(os.path.dirname(__file__), "tcl_parity", "re_annotation.tcl")
pytestmark = pytest.mark.skipif(shutil.which("tclsh") is None or shutil.which("bedtools") is None,
                                reason="tclsh or bedtools not installed")
SAMPLE = [p for p in INPUTS if re.search(r"test_(02|07|08|16|18|24|28|30|33)_", p)][:8]


def _gz(path):
    return io.TextIOWrapper(gzip.GzipFile(path, "wb", mtime=0))


def _el(rnd, pos):
    s = max(0, pos + rnd.randint(-30000, 30000))
    return s, s + rnd.randint(50, 4000)


def _annotations(root, input_file, build):
    rnd = random.Random(input_file + build)
    ann = root / "ann" / "Annotations_Human"
    gdir = ann / "Genes" / build
    gdir.mkdir(parents=True)
    rows = _genes_rows(input_file)
    with _gz(gdir / "ncbiRefSeq.txt.gz") as f:
        f.write("".join("\t".join(r) + "\n" for r in rows))
    ens = [[r[2][3:], r[4], r[5], r[3], r[12], r[1].replace("NM_", "ENST"), r[6], r[7], r[9], r[10]]
           for r in rows if not r[2].endswith("_alt") and rnd.random() < 0.8]
    ens.sort(key=lambda r: (r[0], int(r[1])))
    (gdir / "genes.ENSEMBL.sorted.tmp.bed").write_text("".join("\t".join(r) + "\n" for r in ens))
    genes = sorted({r[12] for r in rows}) + ["UNKNOWN1", "UNKNOWN2"]
    pick = lambda: rnd.choice(genes)  # noqa: E731
    bps = [(re.sub("^chr", "", c), p) for c, p in _breakpoints(input_file)][:150]
    reg = ann / "FtIncludedInSV" / "RegulatoryElements"
    for b in ("GRCh37", "GRCh38"):
        (reg / b).mkdir(parents=True)
    # EnhancerAtlas
    with open(reg / "GRCh37" / "x_EP.txt", "w") as f:
        for c, p in bps:
            for _ in range(rnd.randint(0, 3)):
                s, e = _el(rnd, p)
                f.write(f"chr{c}:{s}-{e}_ENSG0001${pick()}$chr{c}${p}$+\t1.2\n")
        f.write("bad line\n")
    # GeneHancer
    types = ["Enhancer", "Promoter/Enhancer", "promoter/enhancer", "Promoter"]
    el, asso, hg19 = ["chr\telement_start\telement_end\tGHid\tis_elite\tregulatory_element_type"], \
        ["GHid\tsymbol\tcombined_score\tis_elite"], ["chromosome\tstart\tend\tGHid"]
    for i, (c, p) in enumerate(bps):
        for j in range(rnd.randint(0, 2)):
            gid = f"GH{i}_{j}"
            s, e = _el(rnd, p)
            el.append(f"{c}\t{s}\t{e}\t{gid}\t1\t{rnd.choice(types)}")
            for _ in range(rnd.randint(0, 3)):
                asso.append(f"{gid}\t{pick()}\t1.0\t0")
            if rnd.random() < 0.8:
                s, e = _el(rnd, p)
                hg19.append(f"chr{c}\t{s}\t{e}\t{gid}")
    (reg / "GeneHancer_AnnotSV_elements_v5.25.txt").write_text("\n".join(el) + "\n")
    (reg / "GeneHancer_AnnotSV_gene_association_scores_v5.25.txt").write_text("\n".join(asso) + "\n")
    (reg / "GeneHancer_AnnotSV_hg19_v5.25.txt").write_text("\n".join(hg19) + "\n")
    # miRTargetLink
    mi = ["seqnames\tstart\tend\tstrand\tID\tAlias\tName\tspecies\ttarget_gene\tconfidence"]
    for c, p in bps:
        if rnd.random() < 0.5:
            s, e = _el(rnd, p)
            for _ in range(rnd.randint(1, 3)):
                mi.append(f"chr{c}\t{s}\t{e}\t-\tM\tM\thsa-x\thsa\t{pick()}\tFunctional MTI")
    (reg / "Validated_miRNA-gene_pairs_hsa_miRBase_v22.1_GRCh38_location_augmented.tsv").write_text("\n".join(mi) + "\n")
    # ABC: raw GRCh37 download + unsorted formatted GRCh38 files
    abc = ["chr\tstart\tend\tname\tclass\tactivity_base\tTargetGene\tTargetGeneTSS"]
    unsorted = []
    for c, p in bps:
        for _ in range(rnd.randint(0, 2)):
            s, e = _el(rnd, p)
            abc.append(f"chr{c}\t{s}\t{e}\tgenic|x\tgenic\t0.4\t{pick()}\t{p}")
            unsorted.append(f"{c}\t{s}\t{e}\tABC_enhancer\t{pick()};{pick()}")
    with _gz(reg / "GRCh37" / "AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz") as f:
        f.write("\n".join(abc) + "\n")
    rnd.shuffle(unsorted)
    for tx in ("RefSeq", "ENSEMBL"):
        (reg / "GRCh38" / f"ABC_{tx}_GRCh38.sorted.bed").write_text("\n".join(unsorted) + "\n")
    # MPRA
    for name, hdr in (("promoters", "Name\tGenomic coordinates (GRCh37)\tGenomic coordinates (GRCh38)\tTranscript"),
                      ("enhancers", "Name\tGenomic coordinates (GRCh37)\tGenomic coordinates (GRCh38)\tCell line")):
        lines = [hdr]
        for c, p in bps[: 40]:
            if rnd.random() < 0.5:
                s, e = _el(rnd, p)
                s2, e2 = _el(rnd, p)
                suffix = rnd.choice(["", " (rs1)", "+enh"])
                lines.append(f"{pick()}{suffix}\tchr{c}:{s:,}-{e:,}\tchr{c}:{s2:,}-{e2:,}\tx")
        (reg / f"MPRA_{name}.tsv").write_text("\n".join(lines) + "\n")
    # ClinGen, OMIM morbid
    (ann / "Gene-based" / "ClinGen").mkdir(parents=True)
    (ann / "Gene-based" / "OMIM").mkdir(parents=True)
    hi_vals = ["3", "2", "Not yet evaluated", "30", "1"]
    (ann / "Gene-based" / "ClinGen" / "x_ClinGenAnnotations.tsv").write_text(
        "".join(f"{g}\t{rnd.choice(hi_vals)}\t{rnd.choice(hi_vals)}\n" for g in genes if rnd.random() < 0.5))
    with _gz(ann / "Gene-based" / "OMIM" / "x_morbid.tsv.gz") as f:
        f.write("genes\tmorbid\n" + "".join(f"{g}\tyes\n" for g in genes if rnd.random() < 0.3))
    cand = root / "cand.txt"
    cand.write_text(" ".join(g for g in genes if rnd.random() < 0.1) + "\n")
    return str(root / "ann"), str(cand)


def _python(ann, out, input_file, build, tx, sel1, sel2, report, cand):
    c = Config({"annotationsDir": ann, "organism": "Human", "genomeBuild": build, "outputDir": out,
                "outputFile": "x.annotated.tsv", "SVinputFile": input_file, "svtBEDcol": "4",
                "samplesidBEDcol": "-1", "txFile": "", "annotationMode": "both", "tx": tx,
                "SVinputInfo": "1", "SVminSize": "50", "bracketedAltMaxSize": "10", "bracketedAltMode": "ALL",
                "includeCI": "1", "missingGTinSamplesid": "1", "candidateSnvIndelFiles": [],
                "promoterSize": "500", "REselect1": sel1, "REselect2": sel2, "REreport": report,
                "candidateGenesFile": cand, "outputColHeader": ["AnnotSV_ID", "Samples_ID"]})
    try:
        sv = read_sv_input(c)
        check_genes_refseq_file(c)
        check_genes_ensembl_file(c)
        for name, fn in (("checkPromoterFile", rem.check_promoter_file), ("checkEAfiles", rem.check_ea_files),
                         ("checkGHfiles", rem.check_gh_files), ("checkMiRTargetLinkFiles", rem.check_mirtargetlink_files),
                         ("checkABCfiles", rem.check_abc_files), ("checkMPRAfiles", rem.check_mpra_files)):
            fn(c)
            flags = sorted(f"{k} {c[k]}" for k in ("ABCann", "EAann", "GHann", "MPRAann", "miRNAann", "promAnn") if k in c)
            with open(os.path.join(out, "flags.txt"), "a") as f:
                f.write(f"{name} {' '.join(flags)}\n")
        g = genes_annotation(c)
        r = rem.regulatory_elements_annotation(c, g.all_genes_overlapped)
        with open(os.path.join(out, "hits.txt"), "w") as f:
            for k in sorted(r.hi):
                f.write(f"HI {k}={r.hi[k]}\n")
            for k in sorted(r.ts):
                f.write(f"TS {k}={r.ts[k]}\n")
            for k in sorted(r.re):
                f.write(f"RE {k}={r.re[k]}\n")
        organize_annotation(c, g, vcf_header=sv.vcf_header, svlen=sv.svlen, closest_genes=ClosestGenes(c),
                            regulatory=r.re, header_file_to_remove=sv.header_file_to_remove)
    except ConfigExit as e:
        with open(os.path.join(out, "EXIT"), "w") as f:
            f.write(f"{e.code}\n")


CASES = [  # genomeBuild tx REselect1 REselect2 REreport
    ("GRCh38", "RefSeq", "1", "1", "0"),
    ("GRCh38", "ENSEMBL", "0", "0", "1"),
    ("GRCh37", "RefSeq", "0", "1", "0"),
]


@pytest.mark.parametrize("case", CASES, ids="-".join)
@pytest.mark.parametrize("input_file", SAMPLE, ids=lambda p: os.path.relpath(p, TESTS))
def test_regulatory_elements_match_tcl(tmp_path, input_file, case):
    build, tx, sel1, sel2, report = case
    results = []
    for side in ("t", "p"):
        root = tmp_path / side
        root.mkdir()
        ann, cand = _annotations(root, input_file, build)
        work, out = root / "in", root / "out"
        work.mkdir(), out.mkdir()
        inp = str(work / os.path.basename(input_file))
        shutil.copy(input_file, inp)
        if side == "t":
            r = subprocess.run(["tclsh", DRIVER, ROOT, ann, str(out), inp, build, tx, sel1, sel2, report, cand],
                               capture_output=True, text=True)
            assert r.returncode == 0 and "Error" not in r.stdout[-300:], (r.stdout + r.stderr)[-3000:]
        else:
            _python(ann, str(out), inp, build, tx, sel1, sel2, report, cand)
        results.append(root)
    out_p = results[1] / "out" / "x.annotated.tsv"
    header = open(out_p).readline() if out_p.exists() else ""
    n_cols = len(header.rstrip("\n").split("\t"))
    t, p = (_collect(r / "out", n_cols) for r in results)
    assert sorted(t) == sorted(p)
    for name in t:
        assert t[name].split("\n") == p[name].split("\n"), name
    t_files = sorted(os.path.relpath(os.path.join(d, f), results[0]) for sub in ("in", "ann")
                     for d, _, fs in os.walk(results[0] / sub) for f in fs)
    p_files = sorted(os.path.relpath(os.path.join(d, f), results[1]) for sub in ("in", "ann")
                     for d, _, fs in os.walk(results[1] / sub) for f in fs)
    assert t_files == p_files
    for rel in t_files:
        assert open(results[0] / rel, "rb").read() == open(results[1] / rel, "rb").read(), rel
    if "HG00096" in input_file:  # the synthetic data does reach RE_gene
        assert "RE " in open(results[1] / "out" / "hits.txt").read()
