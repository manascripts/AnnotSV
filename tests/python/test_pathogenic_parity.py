"""Compare the Python pathogenic SV and SNV/indel annotations with the Tcl procs.

Covers checkPathogenicFiles (ClinVar, ClinGen gene/region, dbVar, OMIM morbid
parsers, merge with an existing sorted BED, the CHM13 column removal quirk),
checkPathoSNVindelFile and the P_*, po_P_* and P_snvindel_* columns of
OrganizeAnnotation. The output lines are compared whole (not cut to the
header), so the extra empty fields Tcl writes when no column is requested are
checked too.

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
from test_regions_parity import SAMPLE, _small_input
from annotsv import pathogenic
from annotsv.closest_genes import ClosestGenes
from annotsv.config import Config, ConfigExit
from annotsv.genes import check_genes_refseq_file, genes_annotation
from annotsv.sv_input import read_sv_input
from annotsv.write import organize_annotation

DRIVER = os.path.join(os.path.dirname(__file__), "tcl_parity", "pathogenic.tcl")
pytestmark = pytest.mark.skipif(shutil.which("tclsh") is None or shutil.which("bedtools") is None,
                                reason="tclsh or bedtools not installed")

ALL_COLS = ([f"P_{t}_{f}" for t in ("gain", "loss", "ins", "inv") for f in ("phen", "hpo", "source", "coord")]
            + [f"po_P_{t}_{f}" for t in ("gain", "loss") for f in ("phen", "hpo", "source", "coord", "percent")]
            + ["P_snvindel_nb", "P_snvindel_phen"])
COLSETS = {
    "all": ALL_COLS,
    "some": ["P_loss_coord", "P_ins_phen", "po_P_gain_percent", "P_snvindel_nb"],  # po: "^P_" quirk
    "none": [],
}

DISEASES = ["Disease_A", "Syndrome,_type_2", "not_provided", "Disease_A|not_provided", "Cancer|Trait_B",
            "Another_disorder"]
REVSTAT = ["criteria_provided,_multiple_submitters,_no_conflicts", "criteria_provided,_single_submitter",
           "reviewed_by_expert_panel", "no_assertion_criteria_provided", "no_assertion_provided"]
CLNSIG = ["Pathogenic", "Pathogenic/Likely_pathogenic", "Likely_pathogenic", "Benign", "Pathogenic|risk_factor"]


def _gzw(path, text):
    with gzip.GzipFile(path, "wb", mtime=0) as f:
        f.write(text.encode())


def _info(rnd, k, clnvc):
    hpos = ",".join(f"Human_Phenotype_Ontology:HP:{rnd.randint(1, 9999999):07d}" for _ in range(rnd.randint(0, 3)))
    disdb = "|".join(x for x in [hpos, "MONDO:MONDO:0012056,MedGen:C1837873", "MedGen:CN517202"]
                     if x and rnd.random() < 0.8) or "."
    dn_key = "CLNDNINCL" if rnd.random() < 0.1 else "CLNDN"
    dn = rnd.choice(DISEASES) if rnd.random() < 0.95 else ""
    return (f"ALLELEID={rnd.randint(1, 99999) if rnd.random() < 0.97 else ''};CLNDISDB={disdb};"
            f"{dn_key}={dn};CLNHGVS=NC_0:g.1del;CLNREVSTAT={rnd.choice(REVSTAT)};"
            f"CLNSIG={rnd.choice(CLNSIG)};CLNVC={clnvc};GENEINFO=G{k}:1")


def _sv_ref_alt(rnd, length, clnvc):
    seq = "".join(rnd.choice("ACGT") for _ in range(length))
    if clnvc in ("Deletion", "Inversion", "Microsatellite"):
        return "N" + seq, "N"
    return "N", "N" + seq


def _sites(rnd, bps):
    """Pathogenic regions around the SV breakpoints: inside, across, or covering them."""
    for chrom, pos in bps:
        chrom = re.sub("^chr", "", chrom)
        for _ in range(rnd.randint(0, 3)):
            kind = rnd.random()
            length = rnd.choice([1, 20, 49, 50, 51, 300, 4000, 60000])
            if kind < 0.4:
                start = pos + rnd.randint(0, 3000)          # likely inside the SV
            elif kind < 0.8:
                start = max(0, pos - rnd.randint(0, length))  # across the breakpoint
            else:
                start = max(0, pos - length - rnd.randint(0, 500))
            yield chrom, start, length


def _annotations(root, input_file, chm13):
    rnd = random.Random("pathogenic" + input_file + chm13)
    ann, _ = _genes_annotations(root, input_file, "RefSeq")
    base = root / "ann" / "Annotations_Human"
    bps = _breakpoints(input_file)
    genes = sorted({r[12] for r in _genes_rows(input_file)})
    d38 = base / "FtIncludedInSV" / "PathogenicSV" / "GRCh38"
    d38.mkdir(parents=True)

    # OMIM phenotypes (Gene-based) used by ClinGen genes and morbid genes
    omim = base / "Gene-based" / "OMIM"
    omim.mkdir(parents=True)
    _gzw(omim / "20240101_OMIM-1-annotations.tsv.gz", "genes\tOMIM_ID\n" + "".join(
        f"{g}\t{100000 + i}\n" for i, g in enumerate(genes)))
    _gzw(omim / "20240101_OMIM-2-annotations.tsv.gz", "genes\tOMIM_phenotype\tOMIM_inheritance\n" + "".join(
        f"{g}\t{rnd.choice(['{Dis A}, 1 (3)', '[Trait] B/Dis C', 'Syndrome D'])}\tAD\n"
        for g in genes if rnd.random() < 0.6))

    # ClinVar SV
    lines = ["##fileformat=VCFv4.1", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    for k, (chrom, start, length) in enumerate(_sites(rnd, bps)):
        clnvc = rnd.choice(["Deletion", "Duplication", "Insertion", "Inversion", "Microsatellite"])
        ref, alt = _sv_ref_alt(rnd, min(length, 6000), clnvc)
        lines.append(f"{chrom}\t{start + 1}\t{k}\t{ref}\t{alt}\t.\t.\t{_info(rnd, k, clnvc)}")
    _gzw(d38 / "clinvar_20240101.vcf.gz", "\n".join(lines) + "\n")

    # ClinGen gene and region curation lists
    tail = "\tHaploinsufficiency Description\tTriplosensitivity Score\tDate Last Evaluated"

    def clingen(first, rows):
        return ("#ClinGen curation\n#" + first + "\tcytoBand\tGenomic Location\tHaploinsufficiency Score"
                + tail + "\n" + "".join("\t".join(r) + "\n" for r in rows))

    def loc(chrom, start, length):
        if rnd.random() < 0.1:
            return "tbd"
        return rnd.choice(["chr{}:{}-{}", " chr{}: {}-{}", "{}:{}-{}"]).format(chrom, start, start + length)

    sites = list(_sites(rnd, bps))
    g_rows = [[rnd.choice(genes), "1", "1p", loc(*s), rnd.choice(["3", "3", "2", "30", "Not yet evaluated"]), "x",
               rnd.choice(["3", "1", "0", "3"]), "2020"] for s in sites[::2]]
    r_rows = [[f"ISCA-{i}", f"Region {i}", "1p", loc(*s), rnd.choice(["3", "1"]), "x", rnd.choice(["3", "0"]), "2020"]
              for i, s in enumerate(sites[1::2])]
    (d38 / "ClinGen_gene_curation_list_GRCh38.tsv").write_text(clingen("Gene Symbol\tGene ID", g_rows))
    (d38 / "ClinGen_region_curation_list_GRCh38.tsv").write_text(clingen("ISCA ID\tISCA Region Name", r_rows))

    # dbVar
    head = ("#chr\touterMost\toutermost_start\toutermost_stop\tvariant_count\tvariant_type\tmethod\tanalysis"
            "\tplatform\tstudy\tvariant\tclinical_assertion\n")
    for kind in ("deletions", "duplications", "insertions"):
        rows = [f"{c}\tx\t{s}\t{s + (0 if rnd.random() < 0.1 else n)}\t1\tcopy number loss\tm\ta\tp\tst"
                f"\tnssv{rnd.randint(1, 999999)}\tPathogenic\n" for c, s, n in _sites(rnd, bps)]
        _gzw(d38 / f"GRCh38.nr_{kind}.pathogenic.tsv.gz", "#comment\n" + head + "".join(rows) + "\n")

    # OMIM morbid genes
    morbid = [g for g in genes if rnd.random() < 0.3] + ["NOT_A_GENE", ""]
    _gzw(d38 / "20240101_morbid.tsv.gz", "\n".join(morbid) + "\n")

    # An existing sorted BED to merge with, and one that is left as it is
    old = [f"{c}\t{s}\t{s + n}\tOld phen\tHP:0000001 HP:0000002\tCLN:old{i}\t{c}:{s}-{s + n}"
           for i, (c, s, n) in enumerate(_sites(rnd, bps)) if n >= 50]
    old.sort(key=lambda l: (l.split("\t")[0], int(l.split("\t")[1])))
    (d38 / "pathogenic_Gain_SV_GRCh38.sorted.bed").write_text("\n".join(old[::2]) + "\n")
    (d38 / "pathogenic_Inv_SV_GRCh38.sorted.bed").write_text("\n".join(old[1::2]) + "\n")

    # CHM13: drives the column removal of checkPathogenicFiles
    dchm = base / "FtIncludedInSV" / "PathogenicSV" / "CHM13"
    dchm.mkdir(parents=True)
    if chm13 == "all":
        for t in ("Loss", "Gain", "Ins", "Inv"):
            (dchm / f"pathogenic_{t}_SV_CHM13.sorted.bed").write_text("1\t10\t100\t\t\tdbVar:x\t1:10-100\n")
    elif chm13 == "dbvar":
        _gzw(dchm / "CHM13.nr_deletions.pathogenic.tsv.gz", head + "1\tx\t5\t500\t1\tl\tm\ta\tp\tst\tnssv1\tP\n")
        _gzw(dchm / "CHM13.nr_insertions.pathogenic.tsv.gz", head + "2\tx\t5\t500\t1\tl\tm\ta\tp\tst\tnssv2\tP\n")

    # ClinVar SNV/indel (GRCh38, and GRCh37 which is also converted)
    lines = ["##fileformat=VCFv4.1", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    for k, (chrom, start, _) in enumerate(_sites(rnd, bps + bps)):
        n = rnd.choice([0, 0, 1, 5, 48, 49, 60])
        ref, alt = rnd.choice([("A", "G" + "T" * n), ("C" + "A" * n, "C")])
        lines.append(f"{chrom}\t{start + 1}\t{k}\t{ref}\t{alt}\t.\t.\t{_info(rnd, k, 'single_nucleotide_variant')}")
    for build in ("GRCh38", "GRCh37"):
        d = base / "FtIncludedInSV" / "PathogenicSNVindel" / build
        d.mkdir(parents=True)
        _gzw(d / "clinvar_20240101.vcf.gz", "\n".join(lines) + "\n")
    return ann


def _python(ann, out, input_file, mode, cols):
    c = Config({"annotationsDir": ann, "organism": "Human", "genomeBuild": "GRCh38", "outputDir": out,
                "outputFile": "x.annotated.tsv", "SVinputFile": input_file, "svtBEDcol": "4",
                "samplesidBEDcol": "-1", "txFile": "", "annotationMode": mode, "tx": "RefSeq",
                "SVinputInfo": "1", "SVminSize": "50", "bracketedAltMaxSize": "10", "bracketedAltMode": "ALL",
                "includeCI": "1", "missingGTinSamplesid": "1", "candidateSnvIndelFiles": [], "metrics": "us",
                "overlap": "100", "extann": [], "bedtools": "bedtools",
                "outputColHeader": ["AnnotSV_ID", "Samples_ID"] + cols})
    try:
        sv = read_sv_input(c)
        check_genes_refseq_file(c)
        pathogenic.check_pathogenic_files(c)
        pathogenic.check_patho_snvindel_file(c)
        with open(os.path.join(out, "colheader.txt"), "w") as f:
            f.write("\n".join(c["outputColHeader"]) + "\n")
        g = genes_annotation(c)
        organize_annotation(c, g, vcf_header=sv.vcf_header, svlen=sv.svlen, closest_genes=ClosestGenes(c),
                            header_file_to_remove=sv.header_file_to_remove,
                            pathogenic=pathogenic.PathogenicAnnotations(c))
    except ConfigExit as e:
        with open(os.path.join(out, "EXIT"), "w") as f:
            f.write(f"{e.code}\n")


def _collect(out_dir):
    res = {}
    for name in sorted(os.listdir(out_dir)):
        path = os.path.join(out_dir, name)
        if name.endswith(".bash") or os.path.isdir(path):
            continue
        text = open(path).read()
        if name.endswith(".unannotated.tsv"):    # reciprocal breakends come in Tcl hash order
            text = "\n".join(sorted(text.split("\n")))
        res[re.sub(r"^[0-9]{8}-[0-9]{6}_[0-9]+_", "", name)] = text
    return res


RANKING = ["AnnotSV_ranking_score", "AnnotSV_ranking_criteria", "ACMG_class"]

CASES = [  # annotationMode columns CHM13-data
    ("both", "all", "all"),
    ("split", "some", "all"),
    ("full", "all", "dbvar"),
    ("both", "none", "none"),
]


@pytest.mark.parametrize("case", CASES, ids="-".join)
@pytest.mark.parametrize("which", ["small"] + SAMPLE, ids=lambda p: p if p == "small" else os.path.relpath(p, TESTS))
def test_pathogenic_match_tcl(tmp_path, which, case):
    mode, colset, chm13 = case
    input_file = _small_input(tmp_path / "small.bed") if which == "small" else which
    cols = COLSETS[colset]
    results = []
    for side in ("t", "p"):
        root = tmp_path / side
        root.mkdir()
        ann = _annotations(root, input_file, chm13)
        work, out = root / "in", root / "out"
        work.mkdir(), out.mkdir()
        inp = str(work / os.path.basename(input_file))
        shutil.copy(input_file, inp)
        if side == "t":
            r = subprocess.run(["tclsh", DRIVER, ROOT, ann, str(out), inp, "GRCh38", mode, " ".join(cols)],
                               capture_output=True, text=True)
            assert r.returncode == 0 and "Error" not in r.stdout[-300:], (r.stdout + r.stderr)[-3000:]
            tsv = out / "x.annotated.tsv"
            if tsv.exists():
                # Not ported yet: the benign stubs. Each line is compared whole (extra empty fields included).
                lines = tsv.read_text().replace("\t@STUB@", "").split("\n")
                tsv.write_text("\n".join(lines))
        else:
            _python(ann, str(out), inp, mode, cols)
        results.append(root)
    t, p = (_collect(r / "out") for r in results)
    assert sorted(t) == sorted(p)
    for name in t:
        assert t[name].split("\n") == p[name].split("\n"), name
    t_files, p_files = _tree(results[0]), _tree(results[1])
    assert sorted(t_files) == sorted(p_files)
    for rel in t_files:
        assert t_files[rel] == p_files[rel], rel
    out_tsv = t.get("x.annotated.tsv", "")
    if colset == "all" and "EXIT" not in t:
        assert "P_snvindel_phen" in out_tsv and "po_P_loss_percent" in out_tsv
