"""Compare the Python benign SV annotations with the Tcl procs.

Covers checkBenignFiles (gnomAD v2.1 / v4.0, DGV, DDD, 1000g with and without
bcftools, dbVar, ClinGen, ClinVar, IMH, CMRI, PacBio CoLoRS, FGR and HPRC
parsers, merge with an existing sorted BED, the CHM13 column removal quirk),
checkOverlappedGenesBenignFiles and the B_* and po_B_* columns of
OrganizeAnnotation (us and fr metrics). The output lines are compared whole,
so the extra empty fields Tcl writes when no column is requested are checked
too. Error cases (a missing FGR file, a benign SV type with no overlapped gene,
no ENSEMBL genes) compare the exit code and the files left on disk.

Runs only where ``tclsh`` and ``bedtools`` are installed. Every download is
synthetic, built around the breakpoints of each input.
"""
import gzip
import os
import random
import re
import shutil
import subprocess

import pytest

from test_genes_parity import TESTS, ROOT, _annotations as _genes_annotations, _breakpoints, _genes_rows
from test_gene_based_parity import _tree
from test_pathogenic_parity import _collect, RANKING
from test_regions_parity import SAMPLE, _small_input
from annotsv import benign, pathogenic
from annotsv.closest_genes import ClosestGenes
from annotsv.config import Config, ConfigExit
from annotsv.genes import check_genes_refseq_file, genes_annotation
from annotsv.sv_input import read_sv_input
from annotsv.write import organize_annotation

DRIVER = os.path.join(os.path.dirname(__file__), "tcl_parity", "benign.tcl")
pytestmark = pytest.mark.skipif(shutil.which("tclsh") is None or shutil.which("bedtools") is None,
                                reason="tclsh or bedtools not installed")

ALL_COLS = ([f"B_{t}_{f}" for t in ("gain", "loss", "ins", "inv") for f in ("source", "coord", "AFmax")]
            + [f"po_B_{t}_{g}_{f}" for t in ("gain", "loss") for g in ("allG", "someG") for f in ("source", "coord")])
COLSETS = {
    "all": ALL_COLS,
    "some": ["B_loss_coord", "B_inv_AFmax", "po_B_gain_someG_coord"],
    "none": [],
}
AFS = ["0.0500", "0.012", "0.0099", "0.001", "0.0005", "0.3", "1e-02", "0.25", "0.75"]


def _gzw(path, text):
    with gzip.GzipFile(path, "wb", mtime=0) as f:
        f.write(text.encode())


def _sites(rnd, bps, n=(0, 3)):
    """Benign regions around the SV breakpoints: inside, across, or covering them."""
    for chrom, pos in bps:
        chrom = re.sub("^chr", "", chrom)
        for _ in range(rnd.randint(*n)):
            kind = rnd.random()
            length = rnd.choice([20, 49, 50, 51, 300, 4000, 60000, 200000])
            if kind < 0.4:
                start = pos + rnd.randint(0, 3000)
            elif kind < 0.8:
                start = max(1, pos - rnd.randint(0, length))
            else:
                start = max(1, pos - length - rnd.randint(0, 500))
            yield chrom, start, length


def _vcf(rows, info_defs=""):
    return ("##fileformat=VCFv4.2\n" + info_defs + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            + "".join("\t".join(map(str, r)) + "\n" for r in rows))


def _annotations(root, input_file, variant):
    rnd = random.Random("benign" + input_file + variant)
    ann, _ = _genes_annotations(root, input_file, "RefSeq")
    base = root / "ann" / "Annotations_Human"
    bps = _breakpoints(input_file)
    rows = _genes_rows(input_file)
    # ENSEMBL genes (checkOverlappedGenesBenignFiles uses both transcripts)
    gdir = base / "Genes" / "GRCh38"
    if variant != "noensembl":
        bed = sorted({(r[2][3:], int(r[4]), int(r[5]), r[3], r[12]) for r in rows if not r[2].endswith("_alt")},
                     key=lambda r: (r[0], r[1]))
        (gdir / "genes.ENSEMBL.sorted.bed").write_text(
            "".join(f"{c}\t{s}\t{e}\t{st}\t{g}\tENST{i}\t{s}\t{e}\t{s},\t{e},\n" for i, (c, s, e, st, g) in enumerate(bed)))
    d = base / "SVincludedInFt" / "BenignSV" / "GRCh38"
    d.mkdir(parents=True)
    sites = list(_sites(rnd, bps))

    def af():
        return rnd.choice(AFS)

    # gnomAD v4.0 (2 chromosomes files)
    defs = ""
    for name in ("chr1", "chrX"):
        lines = []
        for k, (c, s, n) in enumerate(sites[:: 2 if name == "chr1" else 3]):
            info = [f"END={s + n}", f"SVTYPE={rnd.choice(['DEL', 'DUP', 'INS', 'INV', 'BND', 'CPX', 'CNV', 'INS:ME'])}",
                    f"controls_and_biobanks_AN={rnd.choice([999, 1000, 20000])}",
                    f"controls_and_biobanks_AF={af()}", f"controls_and_biobanks_N_HOMALT={rnd.choice([4, 5, 80])}"]
            if rnd.random() < 0.1:
                info[0], info[1] = info[1], info[0]     # lsearch -start: the order matters
            if rnd.random() < 0.05:
                info.pop(rnd.randrange(len(info)))
            lines.append([c, s + 1, f"gnomAD-SV_v3_{name}_{k}", "N", "<DEL>", ".",
                          rnd.choice(["PASS", "PASS", "LOWQUAL"]), ";".join(["X=1"] + info + ["Y=2"])])
        _gzw(d / f"gnomad.v4.0.sv.{name}.vcf.gz", _vcf(lines, defs))

    # DGV
    head = ("variantaccession\tchr\tstart\tend\tvarianttype\tvariantsubtype\treference\tpubmedid\tmethod\tplatform"
            "\tmergedvariants\tsupportingvariants\tmergedorsample\tfrequency\tsamplesize\tobservedgains"
            "\tobservedlosses\tcohortdescription\tgenes\tsamples\n")
    lines = [f"dgv{k}\t{c}\t{s}\t{s + n}\tCNV\tloss\tref\t1\tm\tp\t\t\tM\t\t{rnd.choice(['500', '501', '2000', ''])}"
             f"\t{rnd.choice([0, 1, 3, 50])}\t{rnd.choice([0, 2, 9, 500])}\t\tG\tS\n" for k, (c, s, n) in enumerate(sites[1::3])]
    (d / "GRCh38_hg38_variants_2020-02-25.txt").write_text(head + "".join(lines) + "\n")

    # DDD
    head = ("#population_cnv_id\tchr\tstart\tend\tdeletion_observations\tdeletion_frequency\tdeletion_standard_error"
            "\tduplication_observations\tduplication_frequency\tduplication_standard_error\tobservations\tfrequency"
            "\tstandard_error\ttype\tsample_size\tstudy\n")
    lines = [f"{k}\t{c}\t{s}\t{s + n}\t1\t{rnd.choice(['0', '0.002', '0.5', '1e-3'])}\t0\t1"
             f"\t{rnd.choice(['0', '0.0011', '0.02'])}\t0\t2\t0.1\t0\t-1\t{rnd.choice([499, 500, 9000])}\tDDD\n"
             for k, (c, s, n) in enumerate(sites[2::3])]
    _gzw(d / "population_cnv_grch38.txt.gz", head + "".join(lines))

    # 1000g (multiallelic sites split by bcftools)
    defs = ('##FILTER=<ID=q10,Description="q">\n' if variant != "badfilter" else "")
    defs += ("##contig=<ID=1>\n##contig=<ID=2>\n##contig=<ID=3>\n##contig=<ID=X>\n##contig=<ID=Y>\n"
            '##INFO=<ID=END,Number=1,Type=Integer,Description="e">\n'
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="t">\n'
            '##INFO=<ID=AF,Number=A,Type=Float,Description="a">\n'
            '##INFO=<ID=EAS_AF,Number=A,Type=Float,Description="a">\n'
            '##INFO=<ID=CS,Number=1,Type=String,Description="c">\n')
    lines = []
    for k, (c, s, n) in enumerate(sites[::2]):
        if c not in ("1", "2", "3", "X", "Y"):
            continue
        alt, svtype = rnd.choice([("<DEL>", "DEL"), ("<DUP>", "DUP"), ("<CN0>", "DEL"), ("<CN2>", "DEL"),
                                  ("<INS:ME:ALU>", "ALU"), ("<INS:ME:LINE1>", "LINE1"), ("<INV>", "INV"),
                                  ("<CN0>,<CN2>", "CNV"), ("A" * 60, "INS"), ("A" * 10, "INS")])
        n_alt = alt.count(",") + 1
        afs = ",".join(rnd.choice(["0.0005", "0.002", "0.4"]) for _ in range(n_alt))
        info = f"AF={afs};EAS_AF={','.join('0.003' for _ in range(n_alt))};CS=x;SVTYPE={svtype}"
        if not svtype.isupper() or rnd.random() < 0.8:
            info += f";END={s + n}"
        if svtype in ("ALU", "LINE1"):
            info = info.split(";END")[0]
        lines.append([c, s + 1, f"esv{k}", "A", alt, ".", rnd.choice(["PASS", "PASS", "q10"]), info])
    lines.sort(key=lambda r: (r[0], r[1]))
    _gzw(d / "ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz", _vcf(lines, defs))

    # dbVar common (deletions, duplications; 25 duplications on one site: po_B "..." cut at 20)
    c0, s0 = re.sub("^chr", "", bps[0][0]), bps[0][1]
    _gzw(d / "GRCh38.nr_deletions.common.bed.gz", "".join(f"chr{c}\t{s}\t{s + n}\tnssv{k}\n"
                                                         for k, (c, s, n) in enumerate(sites[1::4])) + "\n")
    _gzw(d / "GRCh38.nr_duplications.common.bed.gz", "".join(f"chr{c0}\t{max(1, s0 - 100 - k)}\t{s0 + 100 + k}\tdup{k}\n"
                                                             for k in range(25)))

    # ClinGen gene and region curation lists
    tail = "\tHaploinsufficiency Description\tTriplosensitivity Score\tDate Last Evaluated"

    def clingen(first, rws):
        return ("#ClinGen curation\n#" + first + "\tcytoBand\tGenomic Location\tHaploinsufficiency Score"
                + tail + "\n" + "".join("\t".join(r) + "\n" for r in rws))

    def loc(c, s, n):
        return rnd.choice(["tbd", "chr{}:{}-{}", " chr{}: {}-{}", "{}:{}-{}", "chr{}:{}-{}"]).format(c, s, s + n)

    (d / "ClinGen_gene_curation_list_GRCh38.tsv").write_text(clingen("Gene Symbol\tGene ID", [
        [f"G{k}", "1", "1p", loc(*s), rnd.choice(["40", "3", "40"]), "x", rnd.choice(["40", "1"]), "2020"]
        for k, s in enumerate(sites[3::4])]))
    (d / "ClinGen_region_curation_list_GRCh38.tsv").write_text(clingen("ISCA ID\tISCA Region Name", [
        [f"ISCA-{k}", "R", "1p", loc(*s), rnd.choice(["40", "0"]), "x", rnd.choice(["40", "3"]), "2020"]
        for k, s in enumerate(sites[2::4])]))

    # ClinVar
    lines = []
    for k, (c, s, n) in enumerate(sites[::3]):
        clnvc = rnd.choice(["Deletion", "Duplication", "Insertion", "Inversion", "Microsatellite"])
        seq = "A" * min(n, 400)
        ref, alt = ("N" + seq, "N") if clnvc in ("Deletion", "Inversion") else ("N", "N" + seq)
        info = (f"ALLELEID={k};CLNDISDB=.;CLNREVSTAT={rnd.choice(['criteria_provided,_single_submitter', 'reviewed_by_expert_panel', 'no_assertion_criteria_provided', 'x'])};"
                f"CLNSIG={rnd.choice(['Benign', 'Benign/Likely_benign', 'Pathogenic'])};CLNVC={clnvc};GENEINFO=G:1")
        lines.append([c, s + 1, k, ref, alt, ".", ".", info])
    _gzw(d / "clinvar_20240101.vcf.gz", _vcf(lines))

    # IMH (bedpe)
    head = ("#CHROM_A\tSTART_A\tEND_A\tCHROM_B\tSTART_B\tEND_B\tID\tQUAL\tSTRAND_A\tSTRAND_B\tTYPE\tFILTER\tNAME_A"
            "\tREF_A\tALT_A\tNAME_B\tREF_B\tALT_B\tINFO_A\tINFO_B\n")
    lines = [f"chr{c}\t{s}\t{s + 10}\tchr{c if rnd.random() < 0.9 else 'Y'}\t{s + n - 10}\t{s + n}\t{k}\t9\t-\t-"
             f"\t{rnd.choice(['DEL', 'DUP', 'INV', 'INS', 'BND', 'CNV'])}\tPASS\ta\tN\t<X>\tb\tN\t<X>"
             f"\tSVTYPE=X;AF={af()};NSAMP=3\tSVTYPE=X\n" for k, (c, s, n) in enumerate(sites[1::2])]
    _gzw(d / "IMH.callset.public.bedpe.gz", head + "".join(lines))

    # CMRI and PacBio CoLoRS
    lines = [[f"chr{c}", s + 1, f"pb{k}", "N", "<X>", ".", "PASS",
              f"END={rnd.choice([s + n, s + 1, s - 5 if s > 10 else s])};SVTYPE={rnd.choice(['DEL', 'DUP', 'INS', 'INV', 'BND'])};"
              f"SVLEN={rnd.choice([-n, n, 10])};SVN={rnd.choice([499, 500, 1000])};SVF={af()}"]
             for k, (c, s, n) in enumerate(sites[::2])]
    (d / "pb_joint_merged.sv.vcf").write_text(_vcf(lines))
    lines = [[f"chr{c}", s + 1, f"cl{k}", "N", "<X>", ".", "PASS",
              f"SVTYPE={rnd.choice(['DEL', 'DUP', 'INS', 'INV'])};END={rnd.choice([s + n, s + 1])};SVLEN={-n};"
              f"AC=3;AN={rnd.choice([499, 2762])};NS=1;AF={rnd.choice(['0.04', '0.05', '0.3'])};AC_Het=1"]
             for k, (c, s, n) in enumerate(sites[1::2])]
    _gzw(d / "CoLoRSdb.GRCh38.v1.0.0.pbsv.jasmine.vcf.gz", _vcf(lines))

    # HPRC
    lines = [[f"chr{c}", s + 1, f"h{k}", "N", "<X>", ".", rnd.choice(["PASS", "NearContigEnd"]),
              ";".join(rnd.sample([f"SVTYPE={rnd.choice(['DEL', 'DUP', 'INS', 'INV', 'BND'])}", f"END={s + n}",
                                   f"AF={rnd.choice(['0.04', '0.05', '0.5'])}", "SVLEN=3"], rnd.choice([3, 4])))]
             for k, (c, s, n) in enumerate(sites[::3])]
    _gzw(d / "hprc.GRCh38.pbsv.vcf.gz", _vcf(lines))

    # FGR: one benign SV per type over the genes of the first breakpoints (overlappedGenes never empty)
    genes = [r for r in rows if not r[2].endswith("_alt")][:6]
    for t in ("Gain", "Loss", "Ins", "Inv"):
        if variant == "fgr3" and t == "Ins":
            continue
        fl = [f"{r[2]}\t{max(0, int(r[4]) - rnd.randint(0, 9000))}\t{int(r[5]) + rnd.randint(0, 9000)}\tFGR"
              f"\t{r[2][3:]}:{r[4]}-{r[5]}\t{rnd.choice(['0.0500', 'NA', '0.005', '1e-1'])}" for r in genes]
        if variant == "noinv" and t == "Inv":
            fl = ["chr22\t10\t20\tFGR\t22:10-20\t0.5"]
        (d / f"FGR_benign_{t}_SV_GRCh38.bed").write_text("\n".join(fl) + "\n")

    if variant == "noinv":  # only FGR: no gene overlaps the benign Inv SV
        for f in d.iterdir():
            if not f.name.startswith("FGR_"):
                f.unlink()

    # An existing sorted BED to merge with, and a stale overlappedGenes file
    old = sorted({(c, s, s + n) for c, s, n in sites[::5]}, key=lambda x: (x[0], x[1]))
    (d / "benign_Loss_SV_GRCh38.sorted.bed").write_text(
        "".join(f"{c}\t{s}\t{e}\told{i}\t{c}:{s}-{e}\t0.0200\n" for i, (c, s, e) in enumerate(old)))
    (d / "overlappedGenes_RefSeq_in_benign_Loss_SV_GRCh38.tsv").write_text("1:1-2\tSTALE\n")

    # GRCh37: gnomAD v2.1
    d37 = base / "SVincludedInFt" / "BenignSV" / "GRCh37"
    d37.mkdir(parents=True)
    _gzw(d37 / "gnomad_v2.1_sv.sites.bed.gz", "#chrom\tstart\tend\tname\tsvtype\tN_HOMALT\tPOPMAX_AF\tFILTER\n" + "".join(
        f"{c}\t{s}\t{s + n}\tgn{k}\t{rnd.choice(['DEL', 'DUP', 'INS:ME:ALU', 'MCNV', 'INV'])}\t{rnd.choice([4, 5])}"
        f"\t{rnd.choice(['0.0009', '0.002', '0.1'])}\t{rnd.choice(['PASS', 'LOW'])}\n" for k, (c, s, n) in enumerate(sites)))

    # CHM13: drives the column removal of checkBenignFiles
    dchm = base / "SVincludedInFt" / "BenignSV" / "CHM13"
    dchm.mkdir(parents=True)
    if variant != "chm13part":
        for t in ("Loss", "Gain", "Ins", "Inv"):
            (dchm / f"benign_{t}_SV_CHM13.sorted.bed").write_text("1\t10\t100\tx\t1:10-100\t0.5\n")
    else:
        _gzw(dchm / "CHM13.nr_deletions.common.bed.gz", "1\t5\t500\n")
    return ann


def _python(ann, out, input_file, mode, cols, metrics, bcftools):
    c = Config({"annotationsDir": ann, "organism": "Human", "genomeBuild": "GRCh38", "outputDir": out,
                "outputFile": "x.annotated.tsv", "SVinputFile": input_file, "svtBEDcol": "4",
                "samplesidBEDcol": "-1", "txFile": "", "annotationMode": mode, "tx": "RefSeq",
                "SVinputInfo": "1", "SVminSize": "50", "bracketedAltMaxSize": "10", "bracketedAltMode": "ALL",
                "includeCI": "1", "missingGTinSamplesid": "1", "candidateSnvIndelFiles": [], "metrics": metrics,
                "overlap": "100", "extann": [], "bedtools": "bedtools", "bcftools": bcftools,
                "benignAF": "0.01", "minTotalNumber": "500",
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
                            pathogenic=pathogenic.PathogenicAnnotations(c), benign=benign.BenignAnnotations(c))
    except ConfigExit as e:
        with open(os.path.join(out, "EXIT"), "w") as f:
            f.write(f"{e.code}\n")


CASES = [  # annotationMode columns metrics bcftools data-variant
    ("both", "all", "us", "bcftools", "std"),
    ("full", "all", "fr", "nobcftools", "chm13part"),
    ("split", "some", "us", "bcftools", "std"),
    ("both", "none", "us", "bcftools", "std"),
]
ERROR_CASES = [
    ("both", "all", "us", "bcftools", "badfilter"),  # bcftools norm fails: no split, a Tcl error on "0.4,0.4"
    ("both", "all", "us", "bcftools", "fgr3"),
    ("both", "all", "us", "bcftools", "noinv"),
    ("both", "all", "us", "bcftools", "noensembl"),
]


def _run(tmp_path, input_file, case):
    mode, colset, metrics, bcftools, variant = case
    if bcftools == "bcftools" and shutil.which("bcftools") is None:
        pytest.skip("bcftools not installed")
    cols = COLSETS[colset]
    results = []
    for side in ("t", "p"):
        root = tmp_path / side
        root.mkdir()
        ann = _annotations(root, input_file, variant)
        work, out = root / "in", root / "out"
        work.mkdir(), out.mkdir()
        inp = str(work / os.path.basename(input_file))
        shutil.copy(input_file, inp)
        if side == "t":
            r = subprocess.run(["tclsh", DRIVER, ROOT, ann, str(out), inp, "GRCh38", mode, " ".join(cols),
                                metrics, bcftools], capture_output=True, text=True)
            assert r.returncode == 0, (r.stdout + r.stderr)[-3000:]
            tcl_err = r.stderr
        else:
            _python(ann, str(out), inp, mode, cols, metrics, bcftools)
        results.append(root)
    # checkOverlappedGenesBenignFiles leaves its "<clock seconds>_..." tmp file when bedtools fails
    t, p = ({re.sub(r"^[0-9]+_overlappedGenes", "overlappedGenes", k): v for k, v in _collect(r / "out").items()}
            for r in results)
    assert sorted(t) == sorted(p), tcl_err[-2000:]
    for name in t:
        assert t[name].split("\n") == p[name].split("\n"), name
    t_files, p_files = _tree(results[0]), _tree(results[1])
    assert sorted(t_files) == sorted(p_files)
    for rel in t_files:
        a, b = t_files[rel], p_files[rel]
        if "overlappedGenes_" in rel:  # written in Tcl hash order (read back into an array)
            a, b = sorted(a.split(b"\n")), sorted(b.split(b"\n"))
        assert a == b, rel
    t["_err"] = tcl_err[-2000:]
    return t


@pytest.mark.parametrize("case", CASES, ids="-".join)
@pytest.mark.parametrize("which", ["small"] + SAMPLE, ids=lambda p: p if p == "small" else os.path.relpath(p, TESTS))
def test_benign_match_tcl(tmp_path, which, case):
    input_file = _small_input(tmp_path / "small.bed") if which == "small" else which
    t = _run(tmp_path, input_file, case)
    out_tsv = t.get("x.annotated.tsv", "")
    assert "EXIT" not in t, t.get("_err", "")
    if case[1] == "all":
        assert "po_B_loss_someG_coord" in out_tsv
        assert "dbVar;dbVar" in out_tsv and "..." in out_tsv  # po_B lists cut at 20


@pytest.mark.parametrize("case", ERROR_CASES, ids="-".join)
def test_benign_errors_match_tcl(tmp_path, case):
    t = _run(tmp_path, SAMPLE[0], case)
    assert "EXIT" in t
