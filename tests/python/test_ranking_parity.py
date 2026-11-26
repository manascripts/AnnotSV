"""Compare the Python SV ranking (ACMG classification) with the Tcl procs.

Two levels:

* ``test_ranking_procs_match_tcl``: SVprepareRanking, SVrankingLoss/Gain/INS/INV
  and achieveSVrankingLoss/Gain on synthetic annotated lines that walk every
  criterion (sections 1 to 5, Tcl string/number comparisons, "," decimals,
  SV redundancy, custom criteria), then the output order.
* ``test_ranking_output_matches_tcl``: OrganizeAnnotation with the ranking on
  (pathogenic SV, benign SV and a gene file with the gene-based columns),
  whole output compared: ranking columns, rounding, ACMG class,
  ``-rankFiltering``, ``-metrics fr`` and the row order.

Runs only where ``tclsh`` and ``bedtools`` are installed. All data is synthetic.
"""
import os
import random
import re
import shutil
import subprocess

import pytest

from test_genes_parity import TESTS, ROOT
from test_pathogenic_parity import ALL_COLS as P_COLS, _annotations as _patho_annotations
from test_benign_parity import ALL_COLS as B_COLS, _annotations, _collect, _tree
from test_regions_parity import SAMPLE, _small_input
from annotsv import benign, pathogenic
from annotsv import gene_based as gb
from annotsv.closest_genes import ClosestGenes
from annotsv.config import Config, ConfigExit
from annotsv.genes import check_genes_refseq_file, genes_annotation
from annotsv.ranking import Ranking, tcl_str
from annotsv.sv_input import read_sv_input
from annotsv.write import organize_annotation

HERE = os.path.join(os.path.dirname(__file__), "tcl_parity")
pytestmark = pytest.mark.skipif(shutil.which("tclsh") is None or shutil.which("bedtools") is None,
                                reason="tclsh or bedtools not installed")

GENE_COLS = ["HI", "TS", "GnomAD_pLI", "LOEUF_bin", "DDD_HI_percent", "OMIM_morbid"]
PHENO_COLS = ["Exomiser_gene_pheno_score", "PhenoGenius_specificity"]
BASE = ["AnnotSV_ID", "SV_chrom", "SV_start", "SV_end", "SV_length", "SV_type", "Annotation_mode", "Gene_name",
        "Gene_count", "Tx", "Overlapped_CDS_length", "Overlapped_CDS_percent", "Frameshift", "Exon_count",
        "Location", "Location2", "RE_gene", "P_gain_coord", "P_loss_coord", "P_snvindel_nb", "po_P_gain_coord",
        "po_P_loss_coord", "B_gain_coord", "B_loss_coord", "po_B_gain_allG_coord", "po_B_gain_someG_coord",
        "po_B_loss_allG_coord", "po_B_loss_someG_coord"] + GENE_COLS

LOCATIONS = ["txStart-txEnd", "exon2-exon3", "intron1-intron1", "intron1-intron3", "exon3-txEnd", "intron2-txEnd",
             "txStart-exon1", "txStart-intron2", "intron10-intron12", "exon1-intron1", "intron3-exon3"]
LOCATIONS2 = ["5'UTR", "5'UTR-CDS", "CDS-3'UTR", "3'UTR", "CDS", "5'UTR-3'UTR", "UTR", "CDS-UTR"]
RE_GENES = ["", "", "G1 (HI=3/EX=0.8/RE=EA)", "G2 (EX=0,6/morbid/RE=GH);G3 (EX=0.75)", "G4 (HI=3/RE=EA)",
            "G5 (EX=0.5);G6 (EX=0.4)", "G7 (TS=3/EX=0.0000)", "G8 (EX=0.9)", "G9 (PG=A)"]
SVTYPES = ["DEL", "DEL", "DUP", "DUP", "<CN0>", "<CN3>", "DEL_ALU", "INS", "INV", "TRA", "CNV", "LOSS", "gain"]


def _value(rnd, choices):
    return rnd.choice(choices)


def _lines(rnd, cols, n_sv):
    """Annotated full/split lines for ``n_sv`` SV, in the order OrganizeAnnotation ranks them."""
    out = []
    for k in range(n_sv):
        svtype = rnd.choice(SVTYPES)
        start = rnd.randint(1, 10_000_000)
        aid = f"{rnd.choice(['1', '2', 'X'])}_{start}_{start + rnd.randint(50, 900_000)}_{svtype}_1"
        nb_split = rnd.choice([0, 1, 1, 2, 3, 5])
        nb_genes = rnd.choice(["0", "1", "2", "24", "25", "34", "35", "36", "49", "50", "120", str(nb_split)])
        # mostly coherent full lines, sometimes a "deletion within a known pathogenic region" etc.
        full = {"Annotation_mode": "full", "Gene_count": nb_genes, "RE_gene": _value(rnd, RE_GENES)}
        for c in ("P_gain_coord", "P_loss_coord", "po_P_gain_coord", "po_P_loss_coord", "B_gain_coord",
                  "B_loss_coord", "po_B_gain_allG_coord", "po_B_gain_someG_coord", "po_B_loss_allG_coord",
                  "po_B_loss_someG_coord"):
            full[c] = "1:10-20" if rnd.random() < 0.2 else ""
        rows = [full]
        if rnd.random() < 0.05:
            rows.append(dict(full))           # SV redundancy in the input file
        for j in range(nb_split):
            rows.append({
                "Annotation_mode": "split", "Gene_name": f"G{rnd.randint(0, 30)}", "Tx": f"NM_{j}",
                "HI": _value(rnd, ["3", "3", "2", "1", "", "40"]), "TS": _value(rnd, ["3", "", "1"]),
                "OMIM_morbid": _value(rnd, ["yes", "", "", "no"]),
                "GnomAD_pLI": _value(rnd, ["0.95", "0,91", "0.9", "0.5", "", "1", "0"]),
                "LOEUF_bin": _value(rnd, ["0", "3", "", "9"]),
                "DDD_HI_percent": _value(rnd, ["5", "10", "10,5", "50", "", "9.99"]),
                "Exomiser_gene_pheno_score": _value(rnd, ["0.75", "0,6", "0.4", "", "-1.0", "0.7", "0.5", "0,7"]),
                "PhenoGenius_specificity": _value(rnd, ["A", "B", "C", "", ""]),
                "Location": _value(rnd, LOCATIONS), "Location2": _value(rnd, LOCATIONS2),
                "Overlapped_CDS_length": _value(rnd, ["0", "30", "31", "", "0"]),
                "Overlapped_CDS_percent": _value(rnd, ["0", "5", "10", "10,5", "50", "", "9.9"]),
                "Frameshift": _value(rnd, ["yes", "no", "NA"]), "Exon_count": _value(rnd, ["3", "2", "12", "1"]),
                "P_snvindel_nb": _value(rnd, ["", "", "2"]),
            })
        if rnd.random() < 0.1 and nb_split:
            rows = rows[1:] + rows[:1]        # split lines before the full line
        for r in rows:
            vals = {"AnnotSV_ID": aid, "SV_chrom": aid.split("_")[0], "SV_start": str(start), "SV_type": svtype}
            vals.update(r)
            out.append("\t".join(vals.get(c, "") for c in cols))
    return out


def _python_procs(header, lines, svinput, overrides):
    c = Config({"organism": "Human", "SVinputFile": svinput, "svtTSVcol": "5", "svtBEDcol": "4", "ranking": "1"})
    c.ranking_criteria.update(overrides)
    r = Ranking(c)
    r.prepare(header)
    out = [f"ranking={c['ranking']} svtBEDcol={c['svtBEDcol']} svtTSVcol={c['svtTSVcol']}"]
    order = []
    for line in lines:
        ls = line.split("\t")
        if ls[0] not in order:
            order.append(ls[0])
        if r.enabled:
            r.rank(line, ls[5])
    if r.enabled:
        r.achieve_all(order)
        for aid in r.sort_ids(order):
            e = r.expl.get(aid, {}).get("", "")
            out.append(f"{aid}\t{tcl_str(r.score[aid])}\t{e}")
    return "\n".join(out) + "\n"


PROC_CASES = [  # seed, columns, SVinputFile, criteria overrides
    ("all", BASE + PHENO_COLS, "in.bed", {}),
    ("noPheno", BASE, "in.bed", {}),
    ("shuffled", None, "in.vcf.gz", {}),
    ("custom", BASE + PHENO_COLS, "in.bed", {"2A_Loss": "1", "1A_Loss": "0", "5F_Loss": "2", "3A_Gain": "-0.5",
                                             "5G_Gain": "0.333", "2E-4_Loss": "0.25"}),
    ("missingHI", [c for c in BASE if c != "HI"] + ["DDD_HI_percent"], "in.bed", {}),
    ("missingDDD", [c for c in BASE if c != "DDD_HI_percent"], "in.bed", {}),
]


@pytest.mark.parametrize("seed", range(6))
@pytest.mark.parametrize("case", PROC_CASES, ids=lambda c: c[0])
def test_ranking_procs_match_tcl(tmp_path, case, seed):
    name, cols, svinput, overrides = case
    rnd = random.Random(f"{name}{seed}")
    if cols is None:
        cols = BASE[:6] + rnd.sample(BASE[6:] + PHENO_COLS + ["user_HI_x", "TS_extra"], len(BASE) - 6 + 4)
    header = "\t".join(cols)
    lines = _lines(rnd, cols, rnd.choice([5, 40, 200]))
    inp = tmp_path / "lines.tsv"
    inp.write_text(header + "\n" + "\n".join(lines) + "\n")
    out = tmp_path / "tcl.txt"
    ov = " ".join(f"{k} {v}" for k, v in overrides.items())
    r = subprocess.run(["tclsh", os.path.join(HERE, "ranking.tcl"), ROOT, str(inp), str(out), svinput, ov],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-3000:]
    expected = out.read_text()
    got = _python_procs(header, lines, svinput, overrides)
    assert got.split("\n") == expected.split("\n")
    if name in ("all", "custom"):
        assert expected.count("\n") > 3


def _gene_file(root, seed):
    rnd = random.Random(f"genes{seed}")
    rows = ["genes\t" + "\t".join(GENE_COLS + PHENO_COLS)]
    for k in range(0, 400):
        for g in (f"G{k}", f"F{k}"):
            if rnd.random() < 0.85:
                rows.append("\t".join([g, rnd.choice(["3", "3", "2", "1", "0", "40", ""]),
                                       rnd.choice(["3", "1", ""]),
                                       rnd.choice(["0.95", "0.9", "0.5", "", "1", "0.123"]),
                                       rnd.choice(["0", "3", "9", ""]),
                                       rnd.choice(["5", "10", "10.5", "50", "", "9.99"]),
                                       rnd.choice(["yes", "", "yes"]),
                                       rnd.choice(["0.75", "0.6", "0.4", "", "0.7", "0.5"]),
                                       rnd.choice(["A", "B", "C", ""])]))
    path = root / "genes_ranking.tsv"
    path.write_text("\n".join(rows) + "\n")
    return str(path)


def _python(ann, out, input_file, mode, cols, metrics, ext, rank_filtering, svt):
    c = Config({"annotationsDir": ann, "organism": "Human", "genomeBuild": "GRCh38", "outputDir": out,
                "outputFile": "x.annotated.tsv", "SVinputFile": input_file, "svtBEDcol": svt,
                "svtTSVcol": "-1" if svt == "-1" else str(int(svt) + 2),
                "samplesidBEDcol": "-1", "txFile": "", "annotationMode": mode, "tx": "RefSeq",
                "SVinputInfo": "1", "SVminSize": "50", "bracketedAltMaxSize": "10", "bracketedAltMode": "ALL",
                "includeCI": "1", "missingGTinSamplesid": "1", "candidateSnvIndelFiles": [], "metrics": metrics,
                "overlap": "100", "extann": [ext], "geneBasedAnn": "1", "bedtools": "bedtools",
                "bcftools": "bcftools", "benignAF": "0.01", "minTotalNumber": "500",
                "rankFiltering": rank_filtering.split(), "ranking": "1",
                "outputColHeader": ["AnnotSV_ID", "Samples_ID"] + cols})
    try:
        sv = read_sv_input(c)
        check_genes_refseq_file(c)
        pathogenic.check_pathogenic_files(c)
        pathogenic.check_patho_snvindel_file(c)
        benign.check_benign_files(c)
        benign.check_overlapped_genes_benign_files(c)
        with open(os.path.join(out, "colheader.txt"), "w") as f:
            f.write("\n".join(c["outputColHeader"]) + "\n")
        g = genes_annotation(c)
        organize_annotation(c, g, vcf_header=sv.vcf_header, svlen=sv.svlen, closest_genes=ClosestGenes(c),
                            header_file_to_remove=sv.header_file_to_remove,
                            pathogenic=pathogenic.PathogenicAnnotations(c), benign=benign.BenignAnnotations(c),
                            gene_based=gb.GeneBasedAnnotations(c))
    except ConfigExit as e:
        with open(os.path.join(out, "EXIT"), "w") as f:
            f.write(f"{e.code}\n")


OUT_CASES = [  # annotationMode metrics rankFiltering svtBEDcol
    ("both", "us", "1 2 3 4 5 NA", "svt"),
    ("full", "fr", "1 2 3 4 5 NA", "svt"),
    ("split", "us", "3 4 5", "svt"),
    ("both", "us", "1 2 NA", "-1"),
]


@pytest.mark.parametrize("case", OUT_CASES, ids="-".join)
@pytest.mark.parametrize("which", ["small"] + SAMPLE, ids=lambda p: p if p == "small" else os.path.relpath(p, TESTS))
def test_ranking_output_matches_tcl(tmp_path, which, case):
    if shutil.which("bcftools") is None:
        pytest.skip("bcftools not installed")
    mode, metrics, rank_filtering, svt = case
    if svt == "-1" and which != "small" and which.endswith(".vcf"):
        pytest.skip("svtBEDcol is not used with a VCF input")
    input_file = _small_input(tmp_path / "small.bed") if which == "small" else which
    if svt != "-1":
        # SV_type column of the formatted BED ("Samples_ID" = NA is inserted at index 3)
        svt = "5" if which == "small" else "4"
    cols = P_COLS + B_COLS + GENE_COLS + PHENO_COLS
    results = []
    for side in ("t", "p"):
        root = tmp_path / side
        root.mkdir()
        ann = _annotations(root, input_file, "std")
        # The pathogenic SV / SNV-indel data, built aside (same genes) and merged in
        pa = root / "pa"
        pa.mkdir()
        pa_ann = _patho_annotations(pa, input_file, "all")
        for sub in ("FtIncludedInSV", "Gene-based"):
            src = os.path.join(pa_ann, "Annotations_Human", sub)
            if os.path.isdir(src):
                shutil.copytree(src, os.path.join(ann, "Annotations_Human", sub), dirs_exist_ok=True)
        shutil.rmtree(pa)
        # GRCh37 SNV/indels are not used here. With the VCF input, no synthetic ClinVar SNV/indel is
        # kept for GRCh38 either, and an empty conversion is a Tcl error in checkPathoSNVindelFile.
        snv_dir = os.path.join(ann, "Annotations_Human", "FtIncludedInSV", "PathogenicSNVindel")
        shutil.rmtree(snv_dir if input_file.endswith(".vcf") else os.path.join(snv_dir, "GRCh37"),
                      ignore_errors=True)
        ext = _gene_file(root, which)
        work, out = root / "in", root / "out"
        work.mkdir(), out.mkdir()
        inp = str(work / os.path.basename(input_file))
        shutil.copy(input_file, inp)
        if side == "t":
            r = subprocess.run(["tclsh", os.path.join(HERE, "ranking_organize.tcl"), ROOT, ann, str(out), inp,
                                "GRCh38", mode, " ".join(cols), metrics, "bcftools", ext, rank_filtering, svt],
                               capture_output=True, text=True)
            assert r.returncode == 0, (r.stdout + r.stderr)[-3000:]
            assert "EXIT" not in os.listdir(out), (r.stdout + r.stderr)[-3000:]
        else:
            _python(ann, str(out), inp, mode, cols, metrics, ext, rank_filtering, svt)
        results.append(root)
    t, p = (_collect(r / "out") for r in results)
    assert sorted(t) == sorted(p)
    for name in t:
        assert t[name].split("\n") == p[name].split("\n"), name
    t_files, p_files = _tree(results[0]), _tree(results[1])
    assert sorted(t_files) == sorted(p_files)
    out_tsv = t["x.annotated.tsv"]
    assert out_tsv.split("\n")[0].endswith("AnnotSV_ranking_score\tAnnotSV_ranking_criteria\tACMG_class")
    if svt != "-1":
        assert re.search(r"full=[1-5]\n|\t1[AB] \(cf Gene_count", out_tsv)  # the SV are ranked
