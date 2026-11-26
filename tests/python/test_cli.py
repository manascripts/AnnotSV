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


def test_bad_option_exit_code(capsys, install_dir, sv_bed):
    assert cli.main(["-SVinputFile", sv_bed, "-foo", "1"], install_dir=install_dir, check_tools=False) == 2
    assert '"foo" option not known.' in capsys.readouterr().out


def test_config_run_lists_arguments(capsys, install_dir, sv_bed):
    assert cli.main(["-SVinputFile", sv_bed], install_dir=install_dir, check_tools=False) == 1
    assert "AnnotSV has been run with these arguments" in capsys.readouterr().out
