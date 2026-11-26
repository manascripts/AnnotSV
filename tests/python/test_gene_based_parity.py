"""Compare the Python gene-based annotations with the Tcl procs.

Covers the check procs that build the Gene-based files (GenCC, NCBI/HGNC,
memorizeGeneNameAlias, OMIM, morbid genes, HI, ClinGen, LOEUF/pLI, ExAC gene
and CNV intolerance), ExternalAnnotations, the gene-based output columns of
full and split lines (us and fr metrics), isMorbid, isCandidate,
fromOMIM*Phenotype and the candidate genes filtering.

Runs only where ``tclsh`` and ``bedtools`` are installed. The downloads are
synthetic, built with the genes of the synthetic genes table of each input.
"""
import gzip
import io
import os
import random
import re
import shutil
import subprocess

import pytest

from test_genes_parity import INPUTS, TESTS, ROOT, _collect, _genes_rows
from annotsv import gene_based as gb
from annotsv.closest_genes import ClosestGenes
from annotsv.config import Config, ConfigExit, collect_external_gene_files
from annotsv.genes import check_genes_ensembl_file, check_genes_refseq_file, genes_annotation
from annotsv.regulatory_elements import _gene_lists
from annotsv.sv_input import read_sv_input
from annotsv.write import organize_annotation

DRIVER = os.path.join(os.path.dirname(__file__), "tcl_parity", "gene_based.tcl")
pytestmark = pytest.mark.skipif(shutil.which("tclsh") is None or shutil.which("bedtools") is None,
                                reason="tclsh or bedtools not installed")
SAMPLE = [p for p in INPUTS if re.search(r"test_(02|07|16|24|30|33)_", p)][:6]
COLS = ["OMIM_ID", "OMIM_phenotype", "OMIM_inheritance", "OMIM_morbid", "OMIM_morbid_candidate",
        "GenCC_disease", "GenCC_moi", "GenCC_classification", "GenCC_pmid", "HGNC_gene_ID", "HI", "TS",
        "DDD_HI_percent", "LOEUF_bin", "pLI_gnomAD", "pLI_ExAC", "synZ_ExAC", "misZ_ExAC", "delZ_ExAC",
        "dupZ_ExAC", "cnvZ_ExAC", "userScore", "userText", "customCol"]


def _gz(path):
    return io.TextIOWrapper(gzip.GzipFile(path, "wb", mtime=0))


def _dec(rnd):
    v = rnd.choice(["0.5", "-1.25", "9.9443e-01", "3.2e+01", "0.10", "1", "NA", "12", "0.000012", ""])
    return v


def _annotations(root, input_file, build):
    rnd = random.Random("gb" + input_file + build)
    ann = root / "ann" / "Annotations_Human"
    gdir = ann / "Genes" / build
    gdir.mkdir(parents=True)
    rows = _genes_rows(input_file)
    # alias genes: a shifted copy of some transcripts under another name
    genes = sorted({r[12] for r in rows if r[12] != "ALT"})
    alias = {g: f"{g}_AL" for g in genes if rnd.random() < 0.3}
    for r in list(rows):
        if r[12] in alias and not r[2].endswith("_alt"):
            a = list(r)
            a[12], a[1] = alias[r[12]], r[1] + "a"
            a[4], a[5] = str(int(r[4]) + 10), str(int(r[5]) + 10)
            rows.append(a)
    with _gz(gdir / "ncbiRefSeq.txt.gz") as f:
        f.write("".join("\t".join(r) + "\n" for r in rows))
    ens = [[r[2][3:], r[4], r[5], r[3], r[12], r[1].replace("NM_", "ENST"), r[6], r[7], r[9], r[10]]
           for r in rows if not r[2].endswith("_alt") and rnd.random() < 0.8]
    ens.sort(key=lambda r: (r[0], int(r[1])))
    (gdir / "genes.ENSEMBL.sorted.tmp.bed").write_text("".join("\t".join(r) + "\n" for r in ens))
    coord = {}
    for r in rows:
        coord.setdefault(r[12], (r[2], int(r[4]), int(r[5])))
    gene_based = ann / "Gene-based"
    pick = lambda: rnd.choice(genes)  # noqa: E731

    # NCBI / HGNC
    d = gene_based / "NCBIandHGNCgeneID"
    d.mkdir(parents=True)
    res = ["HGNC ID\tApproved symbol\tPrevious symbols\tAlias symbols\tNCBI Gene ID"]
    for i, g in enumerate(genes):
        prev = ", ".join(x for x in (alias.get(g, ""), f"OLD{i}" if rnd.random() < 0.3 else "") if x)
        ncbi = "" if rnd.random() < 0.05 else str(1000 + i)
        res.append(f"HGNC:{i}\t{g}\t{prev}\tAKA{i}\t{ncbi}")
    res.append("HGNC ID\tApproved symbol\tPrevious symbols\tAlias symbols\tNCBI Gene ID")
    (d / "results.txt").write_text("\n".join(res) + "\n")

    # OMIM genemap2 + morbidmap
    d = gene_based / "OMIM"
    d.mkdir(parents=True)
    hdr = ("# Chromosome\tGenomic Position Start\tGenomic Position End\tCyto Location\tComputed Cyto Location\t"
           "MIM Number\tGene/Locus And Other Related Symbols\tGene Name\tApproved Gene Symbol\tEntrez Gene ID\t"
           "Ensembl Gene ID\tComments\tPhenotypes\tMouse Gene Symbol/ID")
    phenos = ["Disease A, 618089 (3), Autosomal dominant", "{Susceptibility B}, 1234 (3)",
              "?Thing C, 55 (2), X-linked recessive, Autosomal recessive", "Syndrome D (4), Y-linked",
              "no number here", "Disorder E, 777 (3), X-linked; [Blood group], 12 (3)",
              "Disease/Trait F, 9 (3), Autosomal recessive"]
    gm = ["# Copyright", hdr]
    for k in range(60):
        g = pick()
        c, s, e = coord[g]
        if rnd.random() < 0.15:
            c = "chr99"
        syms = [g] + ([alias[g]] if g in alias and rnd.random() < 0.5 else []) + (["PRMT9"] if rnd.random() < 0.2 else [])
        approved = g if rnd.random() < 0.7 else ""
        ph = "; ".join(rnd.sample(phenos, rnd.randint(1, 3)))
        gm.append(f"{c}\t{s - 5}\t{e}\t1p\t1p\t{100000 + k % 8}\t{', '.join(syms)}\tname\t{approved}\t1\tENSG\t\t{ph}\tm")
    (d / "genemap2.txt").write_text("\n".join(gm) + "\n")
    mm = ["# comment", "# Phenotype\tGene/Locus And Other Related Symbols\tMIM Number\tCyto Location"]
    mphen = ["Disease, 123 (3)", "[Nondisease], 1 (3)", "Foo (1)", "{Susc}, 2 (3)", "?Prov, 3 (3)", "Bar, 4 (2)",
             "Contiguous (4)"]
    for _ in range(50):
        gs = [pick() for _ in range(rnd.randint(1, 3))] + (["UNKNOWN1"] if rnd.random() < 0.2 else [])
        mm.append(f"{rnd.choice(mphen)}\t{', '.join(gs)}\t5\t1p")
    (d / "morbidmap.txt").write_text("\n".join(mm) + "\n")

    # DDD HI
    d = gene_based / "DDD"
    d.mkdir()
    with _gz(d / "HI_Predictions_Version3.bed.gz") as f:
        f.write("track name=x\n")
        for g in genes:
            if rnd.random() < 0.6:
                f.write(f"chr1\t1\t2\t{g}|0.0004|{rnd.randint(0, 100)}.{rnd.randint(0, 99)}%\t0.1\t.\t1\t2\t0,0,0\n")

    # ClinGen
    d = gene_based / "ClinGen"
    d.mkdir()
    cl = ["#ClinGen Gene Curation Results", "#Gene Symbol\tGene ID\tcytoBand\tGenomic Location\t"
          "Haploinsufficiency Score\tHaploinsufficiency Description\tTriplosensitivity Score\tTD"]
    vals = ["3", "2", "1", "0", "30", "40", "Not yet evaluated"]
    for g in genes:
        for _ in range(1 if rnd.random() < 0.9 else 2):
            if rnd.random() < 0.6:
                cl.append(f"{g}\t1\tp\tloc\t{rnd.choice(vals)}\tx\t{rnd.choice(vals)}\ty")
    (d / "ClinGen_gene_curation_list_GRCh38.tsv").write_text("\n".join(cl) + "\n")

    # gnomAD LOEUF / pLI
    d = gene_based / "gnomAD"
    d.mkdir()
    lo = ["gene\ttranscript\toe_lof_upper_bin\tpLI\texac_pLI"]
    for g in genes:
        for _ in range(rnd.randint(0, 3)):
            lo.append(f"{g}\tT\t{rnd.choice(['0', '3', '9', 'NA', '10'])}\t{_dec(rnd) or 'NA'}\t{_dec(rnd) or 'NA'}")
    (d / "gnomad.v2.1.1.lof_metrics.by_gene.txt").write_text("\n".join(lo) + "\n")

    # GenCC
    d = gene_based / "GenCC"
    d.mkdir()
    gc = ['"uuid"\t"gene_curie"\t"gene_symbol"\t"disease_title"\t"moi_title"\t"classification_title"\t'
          '"submitted_as_pmids"']
    mois = ["Autosomal dominant", "Autosomal dominant inheritance with maternal imprinting", "Autosomal recessive",
            "Digenic inheritance", "Mitochondrial", "Semidominant", "Somatic mosaicism", "X-linked recessive",
            "X-linked", "Y-linked inheritance", "Unknown", ""]
    for k in range(80):
        g = pick() if rnd.random() < 0.9 else "NOTAGENE"
        dis = rnd.choice(['"Disease, with comma"', "Plain disease", '""', '"x ""quoted"" y"'])
        pen = " (incomplete penetrance)" if rnd.random() < 0.1 else ""
        pm = rnd.choice(['"123, 456"', "789", "0", '""', '"11 ,22 33"'])
        gc.append(f'"u{k}"\t"HGNC:1"\t"{g}"\t{dis}\t"{rnd.choice(mois)}{pen}"\t'
                  f'"{rnd.choice(["Definitive", "Moderate", "Limited"])}"\t{pm}')
    (d / "submissions-export-tsv").write_text("\n".join(gc) + "\n")

    # ExAC
    d = gene_based / "ExAC"
    d.mkdir()
    ex = ["transcript\tgene\tchr\tsyn_z\tmis_z"]
    for g in genes:
        if rnd.random() < 0.6:
            ex.append(f"T\t{g}\t1\t{_dec(rnd)}\t{_dec(rnd)}")
    (d / "fordist_cleaned_nonpsych_z_pli_rec_null_data.txt").write_text("\n".join(ex) + "\n")
    cnv = ["gene_symbol chr start end del.score dup.score cnv.score flag"]
    for g in genes:
        for _ in range(rnd.randint(0, 2)):
            cnv.append(f"{g} 1 1 2 {rnd.uniform(-3, 3):.4f} {rnd.uniform(-3, 3):.4f} {rnd.uniform(-3, 3):.4f} F")
    (d / "exac-final-cnv.gene.scores071316").write_text("\n".join(cnv) + "\n")

    # A plain .tsv in Gene-based and a user file given with -externalGeneFiles
    d = gene_based / "Custom"
    d.mkdir()
    (d / "custom.tsv").write_text("# comment\nGenes\tcustomCol\n" + "".join(
        f"{g}\t{_dec(rnd)}\n" for g in genes if rnd.random() < 0.5))
    user = ["#gene file", "genes\tuserScore\tuserText"]
    for g in genes:
        if rnd.random() < 0.5:
            user.append(f"{g}\t{_dec(rnd)}\t{rnd.choice(['a{b}', 'plain', '', 'x y'])}")
    user += [f"{genes[0]}\t0.75\tdup", "short", "   ", f"{genes[-1]}\t1\t2\t3"]
    ext = root / "user.tsv"
    ext.write_text("\n".join(user) + "\n")

    cand = root / "cand.txt"
    cand.write_text(" ".join(g for g in genes if rnd.random() < 0.3) + "\n" + "\t".join(rnd.sample(genes, 2)) + "\n")
    return str(root / "ann"), str(cand), str(ext)


def _python(ann, out, input_file, build, tx, mode, metrics, cand_filtering, cand, ext):
    c = Config({"annotationsDir": ann, "organism": "Human", "genomeBuild": build, "outputDir": out,
                "outputFile": "x.annotated.tsv", "SVinputFile": input_file, "svtBEDcol": "4",
                "samplesidBEDcol": "-1", "txFile": "", "annotationMode": mode, "tx": tx,
                "SVinputInfo": "1", "SVminSize": "50", "bracketedAltMaxSize": "10", "bracketedAltMode": "ALL",
                "includeCI": "1", "missingGTinSamplesid": "1", "candidateSnvIndelFiles": [],
                "candidateGenesFile": cand, "candidateGenesFiltering": cand_filtering, "metrics": metrics,
                "hpo": "HP:1", "extann": [], "externalGeneFiles": [ext],
                "outputColHeader": ["AnnotSV_ID", "Samples_ID"] + COLS})
    try:
        sv = read_sv_input(c)
        check_genes_refseq_file(c)
        check_genes_ensembl_file(c)
        gb.check_gencc_gene_file(c, _gene_lists(c))
        gb.check_ncbi_and_hgnc(c)
        gb.memorize_gene_name_alias(c)
        with open(os.path.join(out, "alias.txt"), "w") as f:
            alias = gb._state(c).get("alias", {})
            for k in sorted(alias):
                f.write(f"{k}={' '.join(alias[k])}\n")
        gb.check_omim_file(c)
        gb.check_morbid_file(c)
        gb.check_hi_file(c)
        gb.check_clingen_file(c)
        gb.check_loeuf_file(c)
        gb.check_gene_intolerance_files(c)
        collect_external_gene_files(c)
        g = genes_annotation(c)
        with open(os.path.join(out, "probe.txt"), "w") as f:
            f.write(f"hpo={c['hpo']}\n")
            for k in sorted(g.genes_by_sv):
                for gene in g.genes_by_sv[k].split(";"):
                    f.write(f"{gene} morbid={int(gb.is_morbid(c, gene))} candidate={int(gb.is_candidate(c, gene))} "
                            f"pheno={gb.from_omim_gene_to_phenotype(c, gene)}\n")
            for i in ("100000", "100001", "100002", "100003", "100004", "100005", "999"):
                f.write(f"{i} {gb.from_omim_to_phenotype(c, i)}\n")
        organize_annotation(c, g, vcf_header=sv.vcf_header, svlen=sv.svlen, closest_genes=ClosestGenes(c),
                            header_file_to_remove=sv.header_file_to_remove,
                            gene_based=gb.GeneBasedAnnotations(c))
    except ConfigExit as e:
        with open(os.path.join(out, "EXIT"), "w") as f:
            f.write(f"{e.code}\n")


def _tree(root):
    files = {}
    for sub in ("in", "ann"):
        for d, _, fs in os.walk(root / sub):
            for name in fs:
                p = os.path.join(d, name)
                rel = os.path.relpath(p, root)
                data = gzip.open(p, "rb").read() if p.endswith(".gz") else open(p, "rb").read()
                if "CNV-Zscore" in name:  # written in Tcl hash order (lookup-only cache, see PARITY.md)
                    data = b"\n".join(sorted(data.split(b"\n")))
                files[rel] = data
    return files


CASES = [  # genomeBuild tx annotationMode metrics candidateGenesFiltering
    ("GRCh38", "RefSeq", "both", "us", "0"),
    ("GRCh38", "ENSEMBL", "both", "fr", "1"),
    ("GRCh37", "RefSeq", "full", "fr", "0"),
]


@pytest.mark.parametrize("case", CASES, ids="-".join)
@pytest.mark.parametrize("input_file", SAMPLE, ids=lambda p: os.path.relpath(p, TESTS))
def test_gene_based_match_tcl(tmp_path, input_file, case):
    build, tx, mode, metrics, cand_filtering = case
    results = []
    for side in ("t", "p"):
        root = tmp_path / side
        root.mkdir()
        ann, cand, ext = _annotations(root, input_file, build)
        work, out = root / "in", root / "out"
        work.mkdir(), out.mkdir()
        inp = str(work / os.path.basename(input_file))
        shutil.copy(input_file, inp)
        if side == "t":
            r = subprocess.run(["tclsh", DRIVER, ROOT, ann, str(out), inp, build, tx, mode, metrics, cand_filtering,
                                cand, ext, " ".join(COLS)], capture_output=True, text=True)
            assert r.returncode == 0 and "Error" not in r.stdout[-300:], (r.stdout + r.stderr)[-3000:]
            tsv = out / "x.annotated.tsv"
            if tsv.exists():
                tsv.write_text(re.sub("\t@STUB@", "", tsv.read_text()))
        else:
            _python(ann, str(out), inp, build, tx, mode, metrics, cand_filtering, cand, ext)
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
    for rel in t_files:
        assert t_files[rel] == p_files[rel], rel
    assert "OMIM_ID" in header and "DDD_HI_percent" in header
