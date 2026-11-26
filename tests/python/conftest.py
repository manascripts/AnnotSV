import os

import pytest


@pytest.fixture
def install_dir(tmp_path):
    """A minimal install tree with an empty-ish configfile."""
    etc = tmp_path / "install" / "etc" / "AnnotSV"
    etc.mkdir(parents=True)
    (etc / "configfile").write_text(
        "# comment\n"
        '-annotationMode:   "both"\n'
        "-benignAF:         0.01\n"
        '-hpo:              ""\n'
        "#------------------------\n"
        "# AnnotSV Output columns:\n"
        "#------------------------\n"
        "AnnotSV_ID                      *\n"
        "SV_chrom                        *\n"
        "CytoBand\t\t\t\n"
        "#Closest_left\n"
        "# AnnotSV ranking criteria:\n"
        "1A_Loss 0.5\n"
    )
    return str(tmp_path / "install")


@pytest.fixture
def sv_bed(tmp_path):
    p = tmp_path / "input" / "sv.bed"
    p.parent.mkdir()
    p.write_text("#chrom\tstart\tend\tsvtype\n1\t100\t2000\tDEL\n")
    return str(p)


@pytest.fixture(autouse=True)
def _chdir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("ANNOTSV", raising=False)
