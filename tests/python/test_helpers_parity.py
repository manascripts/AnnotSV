"""Compare the Python ports of the last bash and Tcl scripts with the originals.

* ``searchForAFreePortNumber.sh``   -> ``helpers.search_for_a_free_port`` (needs ``ss``)
* ``checkPhenoGeniusCliInstall.sh`` -> ``helpers.check_phenogenius_cli_install``
  (fake ``poetry`` and ``git`` on the PATH)
* ``lift_Over_a_BED_file.tcl``      -> ``liftover.main`` (fake ``liftOver``; needs ``tclsh``)

All data is synthetic.
"""
import os
import re
import shutil
import socket
import stat
import subprocess
import sys

import pytest

from test_genes_parity import ROOT
from annotsv import helpers, liftover, tcl

BASH = os.path.join(ROOT, "share", "bash", "AnnotSV")
LIFT_TCL = os.path.join(ROOT, "share", "tcl", "AnnotSV", "Scripts", "lift_Over_a_BED_file.tcl")


def _exe(path, text):
    path.write_text(text)
    path.chmod(path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)


# --------------------------------------------------------------------------
# searchForAFreePortNumber.sh
# --------------------------------------------------------------------------

@pytest.mark.skipif(shutil.which("ss") is None, reason="ss not installed")
def test_free_port_matches_bash():
    with socket.socket() as s:
        s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        try:
            s.bind(("0.0.0.0", 10000))  # the first candidate ("10000" sorts first)
            s.listen(1)
        except OSError:
            pass
        want = tcl.exec_command(["bash", os.path.join(BASH, "searchForAFreePortNumber.sh")])
        got = helpers.search_for_a_free_port()
    assert got == want
    assert got[1] != "10000"


def test_free_port_without_proc(monkeypatch):
    monkeypatch.setattr(helpers, "_ports_in_use", lambda: None)
    failed, port = helpers.search_for_a_free_port()
    assert not failed and 5000 <= int(port) <= 50000


# --------------------------------------------------------------------------
# checkPhenoGeniusCliInstall.sh
# --------------------------------------------------------------------------

FAKE_GIT = """#!/bin/bash
echo "git $*"
if [ "$1" == "clone" ]; then mkdir -p PhenoGeniusCli/.git PhenoGeniusCli/.github; touch PhenoGeniusCli/.gitignore; fi
"""

# FAKE_PG_HELP: "ok" (works), "never" (never works), "afterinstall" (works once installed)
FAKE_POETRY = """#!/bin/bash
echo "poetry $*" >> "$FAKE_LOG"
if [ "$1" == "config" ]; then echo "config out"; echo "config err" >&2; exit 0; fi
if [ "$1" == "install" ]; then echo "installing"; touch installed; exit 0; fi
if [ "$1" == "run" ]; then
  if [ "$FAKE_PG_HELP" == "ok" ] || { [ "$FAKE_PG_HELP" == "afterinstall" ] && [ -e installed ]; }; then
    echo "usage: --hpo_list HPO"
  else
    echo "no such option" >&2; exit 2
  fi
fi
"""


@pytest.mark.parametrize("case", ["nopoetry", "fresh_ok", "installed_ok", "afterinstall", "never"])
def test_phenogenius_install_check_matches_bash(tmp_path, monkeypatch, case):
    results = []
    for side in ("bash", "py"):
        root = tmp_path / side
        tools = root / "tools"
        tools.mkdir(parents=True)
        share = root / "share"
        if case == "nopoetry":
            path = str(tools)
            for name in ("which", "mkdir", "grep", "rm"):
                if shutil.which(name):
                    os.symlink(shutil.which(name), tools / name)
        else:
            _exe(tools / "poetry", FAKE_POETRY)
            path = f"{tools}:{os.environ['PATH']}"
        _exe(tools / "git", FAKE_GIT)
        if case in ("installed_ok", "afterinstall", "never"):
            (share / "python3" / "phenogeniuscli" / "PhenoGeniusCli").mkdir(parents=True)
        log = root / "calls.log"
        env = {"PATH": path, "FAKE_LOG": str(log), "HOME": str(root),
               "FAKE_PG_HELP": {"fresh_ok": "ok", "installed_ok": "ok"}.get(case, case)}
        if side == "bash":
            r = subprocess.run(["/bin/bash", os.path.join(BASH, "checkPhenoGeniusCliInstall.sh"), str(share)],
                               env=env, capture_output=True, text=True)
            out, err = r.stdout.rstrip("\n"), r.stderr.rstrip("\n")
            msg = out + ("\n" + err if out and err else err)
            if not err and r.returncode:
                msg = (msg + "\n" if msg else "") + "child process exited abnormally"
            got = (r.returncode != 0 or err != "", msg)
        else:
            for k, v in env.items():
                monkeypatch.setenv(k, v)
            got = helpers.check_phenogenius_cli_install(str(share))
        files = sorted(os.path.relpath(os.path.join(d, f), share) for d, _, fs in os.walk(share) for f in fs + _)
        logs = {f: (share / f).read_text() for f in files if f.endswith(".log")}
        results.append((got, files, logs, log.read_text() if log.exists() else ""))
    assert results[1] == results[0]


# --------------------------------------------------------------------------
# lift_Over_a_BED_file.tcl
# --------------------------------------------------------------------------

# Shifts the positions by 1000, drops "chrUn" lines (unmapped) and complains on stderr like liftOver
FAKE_LIFTOVER = """#!/usr/bin/env python3
import sys
src, chain, out, unmapped = sys.argv[1:5]
sys.stderr.write("Reading liftover chains\\nMapping coordinates\\n")
with open(src) as f, open(out, "w") as o, open(unmapped, "w") as u:
    for line in f:
        if not line.strip():
            continue
        c = line.rstrip("\\n").split("\\t")
        if c[0] == "chrUn":
            u.write("#Deleted in new\\n" + line)
        else:
            o.write("\\t".join([c[0], str(int(c[1]) + 1000), str(int(c[2]) + 1000)] + c[3:]) + "\\n")
"""

BED = "\n".join([
    "#chrom\tstart\tend\tinfo\tpos",
    "chr2\t500\t900\tB\t2:500-900",
    "chr1\t300\t400\t\tchr1:300-400;x",
    "",
    "chrUn\t10\t20\tlost\t",
    "chr1\t100\t200\tA\tnot:a-position",
    "chrX\t5\t6\t\t\t",
]) + "\n"


def _run_tcl(cwd, args):
    r = subprocess.run(["tclsh", LIFT_TCL] + args, cwd=cwd, capture_output=True, text=True)
    return r.stdout + r.stderr


def _run_py(cwd, args, capsys, monkeypatch):
    monkeypatch.chdir(cwd)
    liftover.main(args)
    c = capsys.readouterr()
    return c.out + c.err


def _norm(text):
    return re.sub(r"\.[0-9]+\.(tmp\.)?bed", ".PID.bed", text)


@pytest.mark.skipif(shutil.which("tclsh") is None, reason="tclsh not installed")
@pytest.mark.parametrize("chain", ["b38-GCA_009914755.4.over.chain", "hg19ToHg38.over.chain"])
def test_liftover_matches_tcl(tmp_path, capsys, monkeypatch, chain):
    tool = tmp_path / "fakeLiftOver"
    _exe(tool, FAKE_LIFTOVER.replace("/usr/bin/env python3", sys.executable))
    res = []
    for side in ("t", "p"):
        d = tmp_path / side
        d.mkdir()
        (d / "in.bed").write_text(BED)
        (d / chain).write_text("")
        args = ["in.bed", "out.bed", chain, str(tool)]
        log = _run_tcl(d, args) if side == "t" else _run_py(d, args, capsys, monkeypatch)
        res.append((_norm(log), (d / "out.bed").read_text(), sorted(os.listdir(d))))
    assert res[1] == res[0]
    assert "chr1\t1100\t1200\tA\tnot:a-position" in res[0][1]


@pytest.mark.skipif(shutil.which("tclsh") is None, reason="tclsh not installed")
def test_liftover_usage_and_tool_error_match_tcl(tmp_path, capsys, monkeypatch):
    tool = tmp_path / "badLiftOver"
    _exe(tool, "#!/bin/bash\necho 'ERROR: bad chain' >&2\nexit 1\n")
    for args in (["in.bed", "out.bed"], ["in.bed", "out.bed", "b38-x.chain", str(tool)]):
        res = []
        for side in ("t", "p"):
            d = tmp_path / f"{side}{len(args)}"
            d.mkdir()
            (d / "in.bed").write_text(BED)
            (d / "out.bed").write_text("old")
            log = _run_tcl(d, args) if side == "t" else _run_py(d, args, capsys, monkeypatch)
            res.append((_norm(log).replace(LIFT_TCL, "AnnotSV-liftover"), sorted(os.listdir(d))))
        assert res[1] == res[0]


def test_liftover_all_mapped(tmp_path, capsys, monkeypatch):
    """The Tcl script stops when nothing is unmapped (grep -c exits 1); the port writes the output."""
    tool = tmp_path / "fakeLiftOver"
    _exe(tool, FAKE_LIFTOVER.replace("/usr/bin/env python3", sys.executable))
    (tmp_path / "in.bed").write_text("chr1\t1\t2\t1:1-2\n")
    _run_py(tmp_path, ["in.bed", "out.bed", "b38-x.chain", str(tool)], capsys, monkeypatch)
    assert (tmp_path / "out.bed").read_text() == "chr1\t1001\t1002\tb38_1:1-2\n"
