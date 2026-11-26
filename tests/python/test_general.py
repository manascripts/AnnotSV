import gzip

from annotsv import general as g
from annotsv.config import Config


def test_file_helpers(tmp_path):
    p = tmp_path / "f.txt"
    g.replace_text_in_file("a\nb", str(p))
    g.write_text_in_file("c", str(p))
    assert p.read_text() == "a\nb\nc\n"
    assert g.content_from_file(str(p)) == "a\nb\nc"
    assert g.lines_from_file(str(p)) == ["a", "b", "c"]
    assert g.first_line_from_file(str(p)) == "a"
    gz = tmp_path / "f.txt.gz"
    with gzip.open(gz, "wt") as f:
        f.write("x\ny\n")
    assert g.lines_from_gz_file(str(gz)) == ["x", "y"]
    assert g.first_line_from_file(str(gz)) == "x"


def test_remove_ref_in_bracketed_alt():
    assert g.remove_ref_in_bracketed_alt("T", "TAAAA[13:5000[") == ("", "AAAA;13:5000;")
    assert g.remove_ref_in_bracketed_alt("T", "]13:5000]AAAAT") == ("", ";13:5000;AAAA")
    assert g.remove_ref_in_bracketed_alt("T", "]22:3000]A") == ("", ";22:3000;")
    assert g.remove_ref_in_bracketed_alt("A", "<DEL>") == ("A", "<DEL>")


def test_annotsv_id_registry():
    reg = g.AnnotSVIDRegistry()
    a = reg.get("1_100_200_DEL", "A", "<DEL>")
    b = reg.get("1_100_200_DEL", "C", "<DEL>")
    assert (a, b) == ("1_100_200_DEL_1", "1_100_200_DEL_2")
    assert reg.get("1_100_200_DEL", "A", "<DEL>") == a
    # BND rescue: REF replaced by N must give the same ID
    x = reg.get("2_5_5_BND", "T", "T[13:5000[")
    assert reg.get("2_5_5_BND", "N", "N[13:5000[") == x


def test_is_not_an_rs():
    assert g.is_not_an_rs("rs123") == 0
    assert g.is_not_an_rs("rs12a") == 1


def test_merge_annotation():
    assert g.merge_annotation(["g1C1\tg1C2"], 2) == "g1C1\tg1C2"
    assert g.merge_annotation(["g1C1\tg1C2", "g2C1\tg2C2", "\t"], 2) == "g1C1;g2C1\tg1C2;g2C2"
    assert g.merge_annotation(["", "x"], 1) == "x"


def test_remove_redundancy():
    assert g.remove_redundancy_without_sorting(["b", "a", "b", "c"]) == ["b", "a", "c"]


def test_normalize_svtype():
    cases = {
        "DEL": "DEL", "<CN0>": "DEL", "loss": "DEL", "DEL_ALU": "DEL",
        "DUP": "DUP", "<CN3>": "DUP", "MCNV": "DUP",
        "<INV>": "INV", "<INS:ME:ALU>": "INS", "LINE1": "INS",
        "TRA": "TRA", "BND": "",
    }
    for raw, expected in cases.items():
        assert g.normalize_svtype(raw) == expected, raw


def test_is_an_empty_file(tmp_path):
    e = tmp_path / "e.vcf"
    e.write_text("##header\n#CHROM\n\n")
    assert g.is_an_empty_file(str(e)) == 1
    v = tmp_path / "v.vcf.gz"
    with gzip.open(v, "wt") as f:
        f.write("#CHROM\n1\t100\n")
    assert g.is_an_empty_file(str(v)) == 0
    assert g.is_an_empty_file(str(tmp_path / "x.txt")) == 1


def test_check_bed_sorts_dedups_and_strips_chr(tmp_path):
    bed = tmp_path / "a.bed"
    bed.write_text("chr2\t50\t60\nchr1\t300\t400\n1\t100\t200\nchr1\t100\t200\nchrUn\t1\t2\n")
    out = g.check_bed(str(bed), output_dir=str(tmp_path))
    assert out == str(tmp_path / "a.formatted.bed")
    assert (tmp_path / "a.formatted.bed").read_text() == "1\t100\t200\n1\t300\t400\n2\t50\t60\n"
    assert not bed.exists()
    assert g.check_bed(str(bed), output_dir=str(tmp_path)) == 1


def test_check_bed_merge_overlap(tmp_path):
    bed = tmp_path / "m.bed"
    bed.write_text("1\t100\t200\tA\n1\t150\t300\tB\n1\t400\t500\tC\n")
    out = g.check_bed(str(bed), output_dir=str(tmp_path), merge_overlap=True)
    assert open(out).read() == "1\t100\t300\tA/B\n1\t400\t500\tC\n"


def test_check_bed_users_header(tmp_path):
    d = tmp_path / "Users" / "GRCh38"
    d.mkdir(parents=True)
    bed = d / "u.bed"
    bed.write_text("#c\ts\te\tname\n1\t1\t2\tx\n")
    g.check_bed(str(bed), output_dir=str(tmp_path), bed_dir=str(d))
    assert (d / "u.header.tsv").read_text() == "#c\ts\te\tname\n"


def _cfg(**kw):
    c = Config(values={"samplesidBEDcol": "-1", "samplesidTSVcol": "-1",
                       "svtBEDcol": "-1", "svtTSVcol": "-1", "ranking": "1"})
    c.values.update(kw)
    return c


def test_check_samplesid_bed_col(tmp_path):
    bed = tmp_path / "s.bed"
    bed.write_text("#h\n1\t1\t2\tDEL\tS1 S2;S3\n")
    c = _cfg(samplesidBEDcol="5")
    g.check_samplesid_bed_col(c, str(bed))
    assert (c["samplesidBEDcol"], c["samplesidTSVcol"]) == ("4", "6")
    assert bed.read_text() == "#h\n1\t1\t2\tDEL\tS1,S2,S3\n"
    c = _cfg(samplesidBEDcol="9")
    g.check_samplesid_bed_col(c, str(bed))
    assert c["samplesidBEDcol"] == "-1"


def test_check_svt_bed_col(tmp_path):
    bed = tmp_path / "s.bed"
    bed.write_text("1\t1\t2\tDEL\n")
    c = _cfg(svtBEDcol="4")
    g.check_svt_bed_col(c, str(bed))
    assert (c["svtBEDcol"], c["svtTSVcol"]) == ("3", "5")
    c = _cfg(svtBEDcol="2")
    g.check_svt_bed_col(c, str(bed))
    assert (c["svtBEDcol"], c["ranking"]) == ("-1", "0")
