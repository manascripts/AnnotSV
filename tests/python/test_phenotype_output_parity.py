"""Compare the PhenoGenius / Exomiser output columns and the Exomiser sort key with Tcl.

OrganizeAnnotation runs with the ranking and with ``PhenoGeniusCli`` /
``Exomiser`` on. The tool results are mocked: the same random
``g_PhenoGeniusCli`` / ``g_Exomiser`` tables (by gene name, and by NCBI gene
ID for PhenoGenius) are given to the Tcl procs and to the Python objects.
The whole output is compared: the 7 columns on the full and split lines, the
Exomiser score used by the ranking (2H, 2I, 5F...) and the row order
(ranking score, then the Exomiser score of the first line of each SV).

Runs only where ``tclsh``, ``bedtools`` and ``bcftools`` are installed. All data is synthetic.
"""
import os
import random
import shutil
import subprocess

import pytest

from test_genes_parity import TESTS, ROOT
from test_pathogenic_parity import ALL_COLS as P_COLS, _annotations as _patho_annotations
from test_benign_parity import ALL_COLS as B_COLS, _annotations, _collect, _tree
from test_regions_parity import SAMPLE, _small_input
from test_ranking_parity import GENE_COLS
import test_ranking_parity
from annotsv import write
from annotsv.exomiser import Exomiser
from annotsv.phenogenius import PhenoGenius

HERE = os.path.join(os.path.dirname(__file__), "tcl_parity")
pytestmark = pytest.mark.skipif(any(shutil.which(t) is None for t in ("tclsh", "bedtools", "bcftools")),
                                reason="tclsh, bedtools or bcftools not installed")

PHENO_OUT_COLS = ["PhenoGenius_score", "PhenoGenius_phenotype", "PhenoGenius_specificity",
                  "Exomiser_gene_pheno_score", "Human_pheno_evidence", "Mouse_pheno_evidence", "Fish_pheno_evidence"]


def _gene_file(root):
    rnd = random.Random("genes-pheno")
    rows = ["genes\t" + "\t".join(GENE_COLS)]
    for k in range(400):
        for g in (f"G{k}", f"F{k}"):
            if rnd.random() < 0.85:
                rows.append("\t".join([g, rnd.choice(["3", "2", "1", "0", ""]), rnd.choice(["3", "1", ""]),
                                       rnd.choice(["0.95", "0.5", "", "1"]), rnd.choice(["0", "3", "9", ""]),
                                       rnd.choice(["5", "10", "50", ""]), rnd.choice(["yes", ""])]))
    path = root / "genes_ranking.tsv"
    path.write_text("\n".join(rows) + "\n")
    return str(path)


def _mock(root, ann, seed):
    """Mock tool results (PG / EX lines) and the NCBI gene IDs of some genes."""
    rnd = random.Random(f"mock{seed}")
    lines, ncbi = [], ["HGNC ID\tApproved symbol\tPrevious symbols\tAlias symbols\tNCBI Gene ID"]
    for k in range(400):
        for g in (f"G{k}", f"F{k}"):
            if rnd.random() < 0.3:
                ncbi.append(f"HGNC:{len(ncbi)}\t{g}\t\t\t{1000 + len(ncbi)}")
                if rnd.random() < 0.5:
                    lines.append(f"PG\t{1000 + len(ncbi) - 1}\t{rnd.choice(['2.5', '0.3'])}\tby id\t"
                                 f"{rnd.choice('ABCD')}")
            if rnd.random() < 0.4:
                lines.append(f"PG\t{g}\t{rnd.choice(['3.2', '0.1', '1'])}\tRenal cysts, Tumor\t"
                             f"{rnd.choice(['A', 'B', 'C', 'D', ''])}")
            if rnd.random() < 0.5:
                score = rnd.choice(["0.8791", "0.7000", "0.6999", "0.5000", "0.4000", "0.0000", "1.0000", "-1.0"])
                lines.append(f"EX\t{g}\t{score}\t{rnd.choice(['', 'Disease A;Brachydactyly'])}\t"
                             f"{rnd.choice(['', 'short digits'])}\t{rnd.choice(['', 'fin'])}")
    d = os.path.join(ann, "Annotations_Human", "Gene-based", "NCBIandHGNCgeneID")
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "results.txt"), "w") as f:
        f.write("\n".join(ncbi) + "\n")
    path = root / "mock.tsv"
    path.write_text("\n".join(lines) + "\n")
    return str(path), lines


def _python(ann, out, inp, mode, cols, metrics, ext, rank_filtering, svt, lines, pg_on, exo_on, monkeypatch):
    """test_ranking_parity._python, with the mocked PhenoGenius / Exomiser objects."""
    objs = {}

    def organize(c, *args, **kwargs):
        c["PhenoGeniusCli"], c["Exomiser"] = pg_on, exo_on
        p, e = PhenoGenius(c), Exomiser(c)
        for L in lines:
            ls = L.split("\t")
            (p if ls[0] == "PG" else e).table[ls[1]] = "\t".join(ls[2:])
        objs["c"] = c
        return orig(c, *args, phenogenius=p, exomiser=e, **kwargs)

    orig = write.organize_annotation
    monkeypatch.setattr(test_ranking_parity, "organize_annotation", organize)
    test_ranking_parity._python(ann, out, inp, mode, cols, metrics, ext, rank_filtering, svt)


CASES = [  # annotationMode metrics rankFiltering PhenoGeniusCli Exomiser
    ("both", "us", "1 2 3 4 5 NA", "1", "1"),
    ("full", "fr", "1 2 3 4 5 NA", "1", "1"),
    ("split", "us", "3 4 5", "1", "1"),
    ("both", "us", "1 2 3 4 5 NA", "0", "1"),
    ("both", "us", "1 2 3 4 5 NA", "1", "0"),
]


@pytest.mark.parametrize("case", CASES, ids="-".join)
@pytest.mark.parametrize("which", ["small"] + SAMPLE, ids=lambda p: p if p == "small" else os.path.relpath(p, TESTS))
def test_phenotype_output_matches_tcl(tmp_path, monkeypatch, which, case):
    mode, metrics, rank_filtering, pg_on, exo_on = case
    input_file = _small_input(tmp_path / "small.bed") if which == "small" else which
    svt = "5" if which == "small" else "4"
    cols = P_COLS + B_COLS + GENE_COLS + PHENO_OUT_COLS
    results = []
    for side in ("t", "p"):
        root = tmp_path / side
        root.mkdir()
        ann = _annotations(root, input_file, "std")
        pa = root / "pa"
        pa.mkdir()
        pa_ann = _patho_annotations(pa, input_file, "all")
        for sub in ("FtIncludedInSV", "Gene-based"):
            src = os.path.join(pa_ann, "Annotations_Human", sub)
            if os.path.isdir(src):
                shutil.copytree(src, os.path.join(ann, "Annotations_Human", sub), dirs_exist_ok=True)
        shutil.rmtree(pa)
        snv_dir = os.path.join(ann, "Annotations_Human", "FtIncludedInSV", "PathogenicSNVindel")
        shutil.rmtree(snv_dir if input_file.endswith(".vcf") else os.path.join(snv_dir, "GRCh37"),
                      ignore_errors=True)
        ext = _gene_file(root)
        mock, lines = _mock(root, ann, os.path.basename(input_file))
        work, out = root / "in", root / "out"
        work.mkdir(), out.mkdir()
        inp = str(work / os.path.basename(input_file))
        shutil.copy(input_file, inp)
        if side == "t":
            r = subprocess.run(["tclsh", os.path.join(HERE, "pheno_organize.tcl"), ROOT, ann, str(out), inp,
                                "GRCh38", mode, " ".join(cols), metrics, "bcftools", ext, rank_filtering, svt,
                                mock, pg_on, exo_on], capture_output=True, text=True)
            assert r.returncode == 0, (r.stdout + r.stderr)[-3000:]
            assert "EXIT" not in os.listdir(out), (r.stdout + r.stderr)[-3000:]
        else:
            _python(ann, str(out), inp, mode, cols, metrics, ext, rank_filtering, svt, lines, pg_on, exo_on,
                    monkeypatch)
        results.append(root)
    t, p = (_collect(r / "out") for r in results)
    assert sorted(t) == sorted(p)
    for name in t:
        assert p[name].split("\n") == t[name].split("\n"), name
    assert sorted(_tree(results[0])) == sorted(_tree(results[1]))
    header = t["x.annotated.tsv"].split("\n")[0]
    assert ("PhenoGenius_score" in header) == (pg_on == "1")
    assert ("Exomiser_gene_pheno_score" in header) == (exo_on == "1")
