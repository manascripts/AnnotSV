"""Port of the bash helpers of ``share/bash/AnnotSV``.

* ``search_for_a_free_port``         -> ``searchForAFreePortNumber.sh``
* ``check_phenogenius_cli_install``  -> ``checkPhenoGeniusCliInstall.sh``

Both return what the Tcl code got from ``catch {exec <script>} Message``, so
the callers keep the Tcl error handling unchanged.
"""

from __future__ import annotations

import os
import shutil
import socket
import subprocess
from typing import Optional, Set, Tuple

PHENOGENIUS_GIT = "https://github.com/kyauy/PhenoGeniusCli.git"
PHENOGENIUS_TAG = "tags/v.1.1.3"  # last PhenoGenius version tested/validated with AnnotSV


# --------------------------------------------------------------------------
# searchForAFreePortNumber.sh
# --------------------------------------------------------------------------

def _ports_in_use() -> Optional[Set[str]]:
    """Local TCP ports of all sockets (``ss -Htan``), None when /proc is not there."""
    ports: Set[str] = set()
    found = False
    for name in ("/proc/net/tcp", "/proc/net/tcp6"):
        try:
            with open(name) as f:
                next(f, None)  # header
                for line in f:
                    fields = line.split()
                    if len(fields) > 1 and ":" in fields[1]:
                        ports.add(str(int(fields[1].rsplit(":", 1)[1], 16)))
            found = True
        except OSError:
            pass
    return ports if found else None


def _can_bind(port: int) -> bool:
    with socket.socket() as s:
        try:
            s.bind(("0.0.0.0", port))
        except OSError:
            return False
    return True


def search_for_a_free_port() -> Tuple[bool, str]:
    """First free TCP port of 5000..50000 in the script's order.

    ``comm -23 <(seq 5000 50000 | sort) <(ports in use | sort -u) | head -n 1``:
    the candidates are sorted as strings ("10000" comes first, "5000" after
    "49999"). Returns (failed, Message) like ``catch {exec ...} Message``.
    Without /proc (e.g. macOS, where ``ss`` does not exist either), a port is
    free when it can be bound.
    """
    in_use = _ports_in_use()
    for port in sorted(str(p) for p in range(5000, 50001)):
        if in_use is not None:
            if port not in in_use:
                return False, port
        elif _can_bind(int(port)):
            return False, port
    return False, ""


# --------------------------------------------------------------------------
# checkPhenoGeniusCliInstall.sh
# --------------------------------------------------------------------------

def _run_to_log(cmd, cwd, log, mode="a") -> int:
    """``cmd >> log 2>&1`` (mode "a") or ``cmd &> log`` (mode "w")."""
    with open(log, mode) as f:
        try:
            return subprocess.run(cmd, cwd=cwd, stdin=subprocess.DEVNULL, stdout=f,
                                  stderr=subprocess.STDOUT).returncode
        except OSError as exc:
            f.write(f"{exc}\n")
            return 127


def _help_has_hpo_list(cwd, log) -> bool:
    """``poetry run python3 phenogenius_cli.py --help &> log`` then ``grep -c -- "--hpo_list" log == 1``."""
    _run_to_log(["poetry", "run", "python3", "phenogenius_cli.py", "--help"], cwd, log, "w")
    with open(log, errors="replace") as f:
        return sum("--hpo_list" in line for line in f) == 1


def check_phenogenius_cli_install(share_dir: str) -> Tuple[bool, str]:
    """Install PhenoGeniusCli under ``<share_dir>/python3/phenogeniuscli`` if needed and check it runs.

    Returns (failed, Message) like ``catch {exec checkPhenoGeniusCliInstall.sh ...} Message``.
    """
    out, err = [], []

    def poetry_config(cwd):
        r = subprocess.run(["poetry", "config", "keyring.enabled", "false"], cwd=cwd,
                           stdin=subprocess.DEVNULL, capture_output=True, text=True)
        out.append(r.stdout)
        err.append(r.stderr)

    def result(code):
        o = "".join(out)
        e = "".join(err)
        o = o[:-1] if o.endswith("\n") else o
        e = e[:-1] if e.endswith("\n") else e
        msg = o
        if e:
            msg = f"{msg}\n{e}" if msg else e
        elif code != 0:
            msg = f"{msg}\nchild process exited abnormally" if msg else "child process exited abnormally"
        return (code != 0 or e != ""), msg

    # Poetry is required
    if shutil.which("poetry") is None:
        out.append("Poetry seems not to be installed\nWARNING: No PhenoGenius annotations available.\n")
        return result(1)

    # PhenoGeniusCli install (if needed)
    base = os.path.join(share_dir, "python3", "phenogeniuscli")
    os.makedirs(base, exist_ok=True)
    src = os.path.join(base, "PhenoGeniusCli")
    if not os.path.isdir(src):
        install_log = os.path.join(base, "PhenoGeniusCli.install.log")
        _run_to_log(["git", "clone", PHENOGENIUS_GIT], base, install_log)
        if os.path.isdir(src):
            _run_to_log(["git", "checkout", PHENOGENIUS_TAG], src, install_log)
            for name in (".git", ".gitattributes", ".github", ".gitignore"):
                p = os.path.join(src, name)
                if os.path.isdir(p):
                    shutil.rmtree(p, ignore_errors=True)
                elif os.path.lexists(p):
                    os.remove(p)
            poetry_config(src)
            _run_to_log(["poetry", "install"], src, os.path.join(base, "poetry_install.log"), "w")
    if not os.path.isdir(src):  # the bash "cd" fails
        err.append(f"cd: {src}: No such file or directory\n")
        return result(1)

    if _help_has_hpo_list(src, os.path.join(base, "PhenoGenius.run.test1.log")):
        return result(0)
    poetry_config(src)
    _run_to_log(["poetry", "install"], src, os.path.join(base, "poetry_install.log"), "w")  # can be run from different servers
    if _help_has_hpo_list(src, os.path.join(base, "PhenoGenius.run.test2.log")):
        return result(0)
    out.append("phenogenius_cli.py not installed\nWARNING: No PhenoGenius annotations available.\n")
    return result(1)
