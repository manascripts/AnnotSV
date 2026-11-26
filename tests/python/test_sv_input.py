import os

from annotsv.config import Config
from annotsv.filtered_vcf import FilteredVCF
from annotsv.sv_input import read_sv_input


def bed_config(tmp_path, bed_text, **kw):
    bed = tmp_path / "in" / "sv.bed"
    bed.parent.mkdir(exist_ok=True)
    bed.write_text(bed_text)
    out = tmp_path / "out"
    out.mkdir(exist_ok=True)
    values = {"SVinputFile": str(bed), "outputDir": str(out), "svtBEDcol": "-1", "samplesidBEDcol": "-1"}
    values.update(kw)
    return Config(values)


def test_bed_with_header_gets_na_samples_column(tmp_path):
    c = bed_config(tmp_path, "#chrom\tstart\tend\ttype\n1\t100\t2000\tDEL\n\n2\t5\t900\tDUP\n")
    state = read_sv_input(c)
    assert c["bedFile"].endswith("in/sv.NA.bed") and c["NAbedFile"] == c["bedFile"]
    assert open(c["bedFile"]).read() == "1\t100\t2000\tDEL\tNA\n2\t5\t900\tDUP\tNA\n"
    assert state.header_file == f"{tmp_path}/out/sv.NA.header.tsv"
    assert open(state.header_file).read() == "#chrom\tstart\tend\ttype\tSamples_ID\n"
    assert c["samplesidBEDcol"] == "5"
    assert not os.path.exists(f"{tmp_path}/out/sv.header.tsv")


def test_bed_without_header_names_user_columns(tmp_path):
    c = bed_config(tmp_path, "1\t100\t2000\tDEL\tS1\tfoo\n", svtBEDcol="4", samplesidBEDcol="5")
    state = read_sv_input(c)
    assert c["bedFile"] == c["SVinputFile"]
    assert open(state.header_file).read() == "SV_chrom\tSV_start\tSV_end\tSV_type\tSamples_ID\tuser#1\n"
    assert not state.add_na_in_bed


def test_filtered_vcf_from_bed_input_is_disabled(tmp_path):
    c = Config({"SVinputFile": "sv.bed", "candidateSnvIndelFiles": ["a.vcf"], "candidateSnvIndelSamples": ["s"]})
    assert FilteredVCF(c).annotate("FULL") is None
    assert c["candidateSnvIndelFiles"] == [] and c["candidateSnvIndelSamples"] == []


def test_filtered_vcf_lists_htz_positions_in_gene(tmp_path):
    vcf = tmp_path / "snv.vcf"
    vcf.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n"
                   "chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t1/1\n"
                   "chr1\t200\t.\tA\tG\t.\tLowQual\t.\tGT\t0|1\t0/1\n"
                   "chr1\t300\t.\tA\t<DEL>\t.\tPASS\t.\tGT\t0/1\t0/1\n"
                   "chr1\t400\t.\tA\tG\t.\tPASS\t.\tGT\t1/1\t0/1\n")
    c = Config({"SVinputFile": "sv.vcf", "candidateSnvIndelFiles": [str(vcf)], "candidateSnvIndelSamples": ["s1", "s2"],
                "SVminSize": "50", "snvIndelPASS": "1", "formatColNumber": 1})
    fv = FilteredVCF(c)
    assert fv.annotate("FULL") == "\t"
    header = "A\tFORMAT\ts1\ts2\tAnnotation_mode"
    line = ["x", "GT", "0/1", "0/0"] + ["z"] * 13
    assert fv.annotate("1", 50, 500, line, header) == "1_100\t"
    assert fv.annotate("2", 50, 500, line, header) == "\t"
