import datetime
import gzip
import os

import pytest

from annotsv import config as cfg
from annotsv.config import ConfigExit


def run(args, install_dir, **kw):
    kw.setdefault("check_tools", False)
    kw.setdefault("log", lambda *a: None)
    return cfg.configure(args, install_dir=install_dir, **kw)


def test_defaults_configfile_and_arguments(install_dir, sv_bed):
    c = run(["-SVinputFile", sv_bed, "-GENOMEBUILD", "GRCh37", "-overlap", "70"], install_dir)
    assert c["genomeBuild"] == "GRCh37"  # case-insensitive name, canonical key
    assert c["organism"] == "Human"
    assert c["overlap"] == "70"
    assert c["benignAF"] == "0.01"
    assert c["rankFiltering"] == ["1", "2", "3", "4", "5", "NA"]
    assert c["outputFile"] == "sv.annotated.tsv"
    assert c["outputDir"] == os.path.abspath(datetime.date.today().strftime("%Y%m%d") + "_AnnotSV")
    assert c["annotationsDir"] == f"{install_dir}/share/AnnotSV"
    assert c["variantconvertDir"] == f"{install_dir}/share/python3/variantconvert/"
    assert c.ranking_criteria["1A_Loss"] == "0.5"
    assert c.ranking_criteria["1B_Loss"] == -0.60
    cols = c["outputColHeader"]
    assert cols[:3] == ["AnnotSV_ID", "SV_chrom", "CytoBand"]
    assert "Closest_left" not in cols
    assert "Samples_ID" in cols and cols.count("SV_chrom") == 1


def test_configfile_next_to_input_wins(install_dir, sv_bed):
    with open(os.path.join(os.path.dirname(sv_bed), "configfile"), "w") as f:
        f.write('-genomeBuild:    "mm10"\n')
    c = run(["-SVinputFile", sv_bed], install_dir)
    assert c["organism"] == "Mouse"


def test_unknown_option(install_dir, sv_bed):
    with pytest.raises(ConfigExit) as e:
        run(["-SVinputFile", sv_bed, "-nope", "1"], install_dir)
    assert e.value.code == 2 and '"nope" option not known.' in e.value.message


def test_double_dash_accepted(install_dir, sv_bed):
    c = run(["--SVinputFile", sv_bed, "--tx", "ENSEMBL"], install_dir)
    assert c["tx"] == "ENSEMBL"


@pytest.mark.parametrize("opt,val", [
    ("benignAF", "0.5"), ("benignAF", ""), ("bracketedAltMaxSize", "-1"), ("vcf", "2"),
    ("tx", "GENCODE"), ("overlap", "0"), ("overlap", "101"), ("minTotalNumber", "50"),
    ("SVminSize", "-5"), ("annotationMode", "x"), ("variantconvertMode", "x"),
    ("bracketedAltMode", "DEL,XYZ"), ("genomeBuild", "hg19"), ("txFile", "/no/such"),
    ("annotationsDir", "/no/such/dir"),
])
def test_bad_values_exit_2(install_dir, sv_bed, opt, val):
    with pytest.raises(ConfigExit) as e:
        run(["-SVinputFile", sv_bed, f"-{opt}", val], install_dir)
    assert e.value.code == 2


def test_input_file_checks(install_dir, tmp_path):
    with pytest.raises(ConfigExit) as e:
        run([], install_dir)
    assert e.value.code == 2
    with pytest.raises(ConfigExit):
        run(["-SVinputFile", "x.txt"], install_dir)
    with pytest.raises(ConfigExit):
        run(["-SVinputFile", "missing.vcf"], install_dir)
    empty = tmp_path / "empty.vcf"
    empty.write_text("#CHROM\n")
    with pytest.raises(ConfigExit) as e:
        run(["-SVinputFile", str(empty)], install_dir)
    assert e.value.code == 0


def test_metrics_fallback_and_hpo(install_dir, sv_bed, capsys):
    c = run(["-SVinputFile", sv_bed, "-metrics", "uk", "-hpo", "HP:0001156;HP:0001363 bad"], install_dir)
    assert c["metrics"] == "us"
    assert c["hpo"] == "HP:0001156,HP:0001363"
    assert "Bad format for the HPO term: bad. Not used." in capsys.readouterr().out


def test_rank_filtering(install_dir, sv_bed):
    c = run(["-SVinputFile", sv_bed, "-rankFiltering", "4-5,1,NA,7"], install_dir)
    assert c["rankFiltering"] == ["1", "4", "5", "NA"]


def test_output_options(install_dir, sv_bed, tmp_path):
    out = tmp_path / "out" / "deep"
    c = run(["-SVinputFile", sv_bed, "-outputFile", f"{out}/res"], install_dir)
    assert c["outputDir"] == str(out) and c["outputFile"] == "res.tsv"
    assert out.is_dir()
    (out / "res.tsv").write_text("x")
    with pytest.raises(ConfigExit):
        run(["-SVinputFile", sv_bed, "-outputFile", f"{out}/res", "-overwrite", "0"], install_dir)


def test_snvindel_samples(install_dir, sv_bed, tmp_path):
    vcf = tmp_path / "snv.vcf.gz"
    with gzip.open(vcf, "wt") as f:
        f.write("##x\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2\tS1\n")
    c = run(["-SVinputFile", sv_bed, "-snvIndelFiles", str(tmp_path / "snv*.vcf.gz")], install_dir)
    assert c["snvIndelFiles"] == [str(vcf)] and c["snvIndelSamples"] == ["S1", "S2"]
    c = run(["-SVinputFile", sv_bed, "-snvIndelFiles", str(vcf), "-snvIndelSamples", "S2,S9"], install_dir)
    assert c["snvIndelSamples"] == ["S2"]
    c = run(["-SVinputFile", sv_bed, "-snvIndelSamples", "S2"], install_dir)
    assert c["snvIndelSamples"] == []


def test_tool_check(install_dir, sv_bed):
    def old_bedtools(cmd):
        return "bedtools v2.20.1" if cmd[0] == "bedtools" else "bcftools 1.17"
    with pytest.raises(ConfigExit) as e:
        run(["-SVinputFile", sv_bed], install_dir, check_tools=True, tool_runner=old_bedtools)
    assert "2.25" in e.value.message

    def missing(cmd):
        raise FileNotFoundError(cmd[0])
    with pytest.raises(ConfigExit):
        run(["-SVinputFile", sv_bed], install_dir, check_tools=True, tool_runner=missing)
    run(["-SVinputFile", sv_bed], install_dir, check_tools=True,
        tool_runner=lambda cmd: "bedtools v2.31.0")


def test_external_gene_files_and_listing(install_dir, sv_bed):
    ann = os.path.join(install_dir, "share", "AnnotSV", "Annotations_Human", "Gene-based", "OMIM")
    os.makedirs(ann)
    for name in ("omim.tsv", "x_DGV_samplesInStudies.tsv", "y.tsv.gz"):
        open(os.path.join(ann, name), "w").close()
    c = run(["-SVinputFile", sv_bed], install_dir)
    cfg.collect_external_gene_files(c)
    names = sorted(os.path.basename(p) for p in c["extann"])
    assert names == ["omim.tsv", "y.tsv.gz"] and c["geneBasedAnn"] == "1"
    listing = cfg.arguments_listing(c)
    assert "\t-SVinputFile " + sv_bed in listing
    assert "\t-rankFiltering 1 2 3 4 5 NA" in listing
    assert "installDir" not in listing and "outputColHeader" not in listing
