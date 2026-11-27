from annotsv import cli


def test_version(capsys, install_dir):
    assert cli.main(["-version"], install_dir=install_dir) == 0
    assert capsys.readouterr().out.strip() == "AnnotSV 3.5.10"


def test_help_and_no_args(capsys, install_dir):
    with open(f"{install_dir}/commandLineOptions.txt", "w") as f:
        f.write("OPTIONS\n")
    assert cli.main([], install_dir=install_dir) == 0
    out = capsys.readouterr().out
    assert "Arguments are missing" in out and "OPTIONS" in out
    assert cli.main(["-help"], install_dir=install_dir) == 0
    assert cli.main(["--help"], install_dir=install_dir) == 0
    assert cli.main(["-h"], install_dir=install_dir) == 0


def test_help_uses_python_style_names():
    text = (
        "       $ANNOTSV/bin/AnnotSV -SVinputFile 'Path' >& AnnotSV.log &\n"
        "-annotationsDir:                Path of the annotations directory\n"
        "-vcf:                           VCF output (-svtBEDcol needs to be defined too)\n"
        "-help:                          More information\n"
    )
    out = cli.python_style_help(text).split("\n")
    assert out[0] == "       annotsv --sv-input-file 'Path' >& AnnotSV.log &"
    assert out[1] == "--annotations-dir               Path of the annotations directory"
    assert out[2] == "--vcf                           VCF output (--svt-bed-col needs to be defined too)"
    assert out[3].startswith("--help                          More")


def test_version_double_dash(capsys, install_dir):
    assert cli.main(["--version"], install_dir=install_dir) == 0
    assert capsys.readouterr().out.strip() == "AnnotSV 3.5.10"


def test_bad_option_exit_code(capsys, install_dir, sv_bed):
    assert cli.main(["-SVinputFile", sv_bed, "-foo", "1"], install_dir=install_dir, check_tools=False) == 2
    assert '"foo" option not known.' in capsys.readouterr().out


def test_config_run_lists_arguments(capsys, install_dir, sv_bed):
    # No annotation bundle in the test install: the genes file check exits with 2, as in Tcl
    assert cli.main(["-SVinputFile", sv_bed], install_dir=install_dir, check_tools=False) == 2
    out = capsys.readouterr().out
    assert "AnnotSV has been run with these arguments" in out
    assert "Please check your install - Exit with error." in out
