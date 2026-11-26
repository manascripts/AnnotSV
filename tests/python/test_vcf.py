import datetime
import os

import pytest

from annotsv.config import Config, ConfigExit
from annotsv.vcf import dichotomy_search, vcfs_to_bed

HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n"
NOW = datetime.datetime(2026, 1, 2, 3, 4, 5)


def make_config(tmp_path, **kw):
    out = tmp_path / "out"
    out.mkdir(exist_ok=True)
    values = {"outputDir": str(out), "outputFile": "sv.annotated.tsv", "SVminSize": "50",
              "bracketedAltMaxSize": "10", "bracketedAltMode": "ALL", "SVinputInfo": "1",
              "includeCI": "1", "missingGTinSamplesid": "1", "candidateSnvIndelFiles": []}
    values.update(kw)
    return Config(values)


def run(tmp_path, body, **kw):
    vcf = tmp_path / "sv.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n" + HEADER + body)
    config = make_config(tmp_path, **kw)
    res = vcfs_to_bed(config, str(vcf), now=NOW)
    bed = [l for l in open(res.bed_file).read().split("\n") if l]
    unann = open(res.unannotated_file).read() if os.path.exists(res.unannotated_file) else ""
    return res, bed, unann


def test_dichotomy_search():
    assert dichotomy_search(5, ["10", "20"]) == -1
    assert dichotomy_search(10, ["10", "20", "30"]) == 0
    assert dichotomy_search(25, ["10", "20", "30"]) == 1
    assert dichotomy_search(30, ["10", "20", "30"]) == 2


def test_angle_bracketed_del(tmp_path):
    res, bed, unann = run(tmp_path, "chr1\t1000\tdel1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=5000;CIPOS=-10,10;CIEND=-5,5\tGT\t0/1\t0/0\n")
    assert bed == ["1\t989\t5005\tDEL\ts1\tdel1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=5000;CIPOS=-10,10;CIEND=-5,5\tGT\t0/1\t0/0"]
    assert res.svlen == {"1_1000_5005_DEL_1": -4000}
    assert res.vcf_header == "SV_type\tSamples_ID\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2"
    assert os.path.basename(res.bed_file) == f"20260102-030405_{os.getpid()}_AnnotSV_inputSVfile.bed"
    assert unann == ""


def test_svinputinfo_0_and_no_ci(tmp_path):
    _, bed, _ = run(tmp_path, "1\t1000\tdel1\tN\t<DEL>\t.\tPASS\tEND=5000;CIPOS=-10,10\tGT\t./.\t1/1\n",
                    SVinputInfo="0", includeCI="0")
    assert bed == ["1\t999\t5000\tDEL\ts1,s2\tN\t<DEL>\tGT\t./.\t1/1"]


def test_small_variants_and_unknown_chrom_are_unannotated(tmp_path):
    body = ("1\t100\t.\tA\tATT\t.\tPASS\t.\tGT\t0/1\t0/1\n"
            "1\t200\t.\tN\t<DUP>\t.\tPASS\tSVLEN=20\tGT\t0/1\t0/1\n"
            "Un_1\t300\t.\tN\t<INV>\t.\tPASS\tEND=5000\tGT\t0/1\t0/1\n"
            "1\t400\t.\tN\tACGT\t.\tPASS\t.\tGT\t0/1\t0/1\n"
            "2\t1000\t.\tN\t<DEL>\t.\tPASS\tEND=9000\tGT\t0/1\t0/1\n")
    _, bed, unann = run(tmp_path, body)
    assert len(bed) == 1 and bed[0].startswith("2\t999\t9000\tDEL")
    assert unann.splitlines() == [
        "1_100___A_ATT: variantLength (2) < SVminSize (50) (line 3)",
        "1_200_200_DUP_N_<DUP>: variantLength (20) < SVminSize (50) (line 4)",
        'Un_1_300_5000_INV_N_<INV>: chromosome "Un_1" unknown (line 5)',
        "1_400___N_ACGT: variantLength (3) < SVminSize (50) (line 6)",
    ]


def test_square_bracketed_del_pair_kept_once(tmp_path):
    body = ("12\t3000\tdel_a\tT\tT[12:5000[\t.\tPASS\tSVTYPE=BND\tGT\t0/1\t0/0\n"
            "12\t5000\tdel_b\tT\t]12:3000]T\t.\tPASS\tSVTYPE=BND\tGT\t0/1\t0/0\n")
    res, bed, unann = run(tmp_path, body)
    assert len(bed) == 1 and bed[0].startswith("12\t2999\t5000\tDEL\ts1\tdel_a")
    assert unann.strip() == "12_3000_5000_DEL_1: reciprocal breakend (line 4)"
    assert res.svlen == {"12_3000_5000_DEL_1": -2000}


def test_square_bracketed_del_rescued_from_second_breakend(tmp_path):
    _, bed, _ = run(tmp_path, "12\t5000\tdel_b\tT\t]12:3000]T\t.\tPASS\tSVTYPE=BND;CIPOS=-1,1;CIEND=-9,9\tGT\t0/1\t0/0\n")
    # includeCI uses the swapped CIPOS/CIEND; REF/ALT stay as read
    assert bed == ["12\t2990\t5001\tDEL\ts1\tdel_b\tT\t]12:3000]T\t.\tPASS"
                   "\tSVTYPE=BND;CIPOS=-9,9;CIEND=-1,1;BNDrescue\tGT\t0/1\t0/0"]


def test_translocation_writes_both_breakends(tmp_path):
    _, bed, _ = run(tmp_path, "17\t198982\ttra_a\tA\tA]2:321681]\t.\tPASS\tSVTYPE=BND\tGT\t0/1\t0/0\n")
    assert sorted(bed) == sorted([
        "17\t198981\t198982\tTRA\ts1\ttra_a\tA\tA]2:321681]\t.\tPASS\tSVTYPE=BND\tGT\t0/1\t0/0",
        "2\t321680\t321681\tTRA\ts1\ttra_a\tA\tA]2:321681]\t.\tPASS\tSVTYPE=BND;BNDrescue\tGT\t0/1\t0/0",
    ])


def test_bracketed_alt_mode_turns_dup_into_bnd(tmp_path):
    _, bed, _ = run(tmp_path, "2\t3000\tdup_a\tT\t]2:5000]T\t.\tPASS\t.\tGT\t0/1\t0/0\n", bracketedAltMode="DEL")
    assert bed[0].startswith("2\t2999\t3000\tBND\t")


def test_insertions_at_same_position_get_distinct_ids(tmp_path):
    body = ("18\t1963003\ta\t.\t<INS>\t.\tPASS\tSVLEN=1512;SVTYPE=INS;END=1963003\tGT\t0/1\t0/0\n"
            "18\t1963003\tb\t.\t<INS>\t.\tPASS\tSVLEN=51;SVTYPE=INS;END=1963003\tGT\t0/1\t0/0\n")
    res, bed, _ = run(tmp_path, body)
    assert [l.split("\t")[7] for l in bed] == ["<INS>_1512", "<INS>_51"]
    assert res.svlen == {"18_1963003_1963003_INS_1": "1512", "18_1963003_1963003_INS_2": "51"}


def test_missing_header_and_empty_result_exit(tmp_path):
    vcf = tmp_path / "nohead.vcf"
    vcf.write_text("1\t1000\t.\tN\t<DEL>\t.\tPASS\tEND=5000\n")
    with pytest.raises(ConfigExit) as exc:
        vcfs_to_bed(make_config(tmp_path), str(vcf))
    assert exc.value.code == 2
    with pytest.raises(ConfigExit) as exc:
        run(tmp_path, "1\t100\t.\tA\tAT\t.\tPASS\t.\tGT\t0/1\t0/1\n")
    assert exc.value.code == 0


def test_sv_gt_absent_disables_compound_htz(tmp_path):
    _, _, _ = run(tmp_path, "1\t1000\t.\tN\t<DEL>\t.\tPASS\tEND=5000\t.\t.\t.\n",
                  candidateSnvIndelFiles=["x.vcf"], candidateSnvIndelSamples=["s1"])
