"""Port of the annotation installer (``bin/INSTALL_annotations.sh``) and of the
data installation targets of the ``Makefile``.

The code itself is installed with ``pip install .`` (it replaces the
``install-executable``, ``install-tcl-toolbox``, ``install-configfile``... targets).
This command installs what pip cannot: the annotation bundles and the external
tools AnnotSV downloads.

    annotsv-install [--annotations-dir DIR] TARGET [TARGET ...]

Targets (Makefile target in brackets):

* ``human``           Annotations_Human_<v>.tar.gz + the Exomiser phenotype data
                      [install-human-annotation, which also runs install-exomiser]
* ``mouse``           Annotations_Mouse_<v>.tar.gz          [install-mouse-annotation]
* ``exomiser``        <v>_phenotype.zip only                [install-exomiser]
* ``exomiser-jar``    exomiser-rest-prioritiser jar          [install-exomiser-resources]
* ``variantconvert``  git clone + pip install of variantconvert [install-variantconvert]

``annotsv-install --annotations-dir ./AnnotSV_annotations human`` is what
``INSTALL_annotations.sh`` did (for use with ``annotsv --annotations-dir``).
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import shutil
import stat
import subprocess
import sys
import tarfile
import urllib.error
import urllib.request
import zipfile
from typing import List, Optional, Sequence

from .config import default_install_dir

# Same defaults as the Makefile variables
HUMAN_VERSION = "3.5"
EXOMISER_VERSION = "2406"
MOUSE_VERSION = "3.4.2"
EXOMISER_JAR_VERSION = "14.1.0"
VC_VERSION = "2.0.1"

ANNOTATIONS_URL = "https://www.lbgi.fr/~geoffroy/Annotations"
EXOMISER_DATA_URL = "https://data.monarchinitiative.org/exomiser/data"
EXOMISER_JAR_URL = "https://github.com/exomiser/Exomiser/releases/download"
VARIANTCONVERT_GIT = "https://github.com/SamuelNicaise/variantconvert.git"


def _log(msg: str = "") -> None:
    print(msg, flush=True)


def download(url: str, dest: str) -> str:
    """``curl -C - -LO url``: download into ``dest``, resuming a partial file."""
    done = os.path.getsize(dest) if os.path.exists(dest) else 0
    req = urllib.request.Request(url, headers={"User-Agent": "AnnotSV-installer"})
    if done:
        req.add_header("Range", f"bytes={done}-")
    try:
        resp = urllib.request.urlopen(req)
    except urllib.error.HTTPError as exc:
        if exc.code == 416:   # already complete
            return dest
        raise
    with resp:
        # A server that ignores Range sends the whole file again
        mode = "ab" if done and getattr(resp, "status", 200) == 206 else "wb"
        with open(dest, mode) as out:
            shutil.copyfileobj(resp, out, 1 << 20)
    return dest


def _fetch(url: str, workdir: str) -> str:
    """Use an already downloaded archive (``make`` would too), else download it."""
    name = url.rsplit("/", 1)[-1]
    for candidate in (os.path.join(workdir, name), name):
        if os.path.isfile(candidate) and _complete(candidate):
            return candidate
    return download(url, os.path.join(workdir, name))


def _complete(path: str) -> bool:
    try:
        if path.endswith(".zip"):
            with zipfile.ZipFile(path) as z:
                return z.testzip() is None
        with tarfile.open(path) as t:
            for _ in t:
                pass
        return True
    except (OSError, EOFError, tarfile.TarError, zipfile.BadZipFile):
        return False


def _chmod_777(path: str) -> None:
    """``chmod -R 777``."""
    mode = stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO
    os.chmod(path, mode)
    for root, dirs, files in os.walk(path):
        for n in dirs + files:
            p = os.path.join(root, n)
            if not os.path.islink(p):
                os.chmod(p, mode)


def _extract_tar(archive: str, dest: str) -> None:
    with tarfile.open(archive) as t:
        if hasattr(tarfile, "data_filter"):
            t.extractall(dest, filter="data")
        else:  # Python < 3.12 without the backported filters
            root = os.path.realpath(dest)
            for m in t.getmembers():
                target = os.path.realpath(os.path.join(dest, m.name))
                if os.path.commonpath([root, target]) != root:
                    raise RuntimeError(f"Unsafe path in {archive}: {m.name}")
            t.extractall(dest)


def _extract_zip(archive: str, dest: str) -> None:
    with zipfile.ZipFile(archive) as z:
        z.extractall(dest)   # zipfile already strips ".." and absolute paths


def install_organism(annotations_dir: str, organism: str, version: str, keep: bool = False) -> None:
    """install-human-annotation / install-mouse-annotation."""
    _log(f"\nInstallation of {organism.lower()} annotation:\n")
    os.makedirs(annotations_dir, exist_ok=True)
    archive = _fetch(f"{ANNOTATIONS_URL}/Annotations_{organism}_{version}.tar.gz", annotations_dir)
    _extract_tar(archive, annotations_dir)
    if not keep:
        os.remove(archive)
    for d in glob.glob(os.path.join(annotations_dir, "Annotations_*")):
        _chmod_777(d)
    _log(f"\n--> {organism} annotation installed")


def install_exomiser(annotations_dir: str, version: str, keep: bool = False) -> None:
    """install-exomiser: Exomiser phenotype data in Annotations_Exomiser/<version>/."""
    _log("\nInstallation of Exomiser data:\n")
    dest = os.path.join(annotations_dir, "Annotations_Exomiser", version)
    os.makedirs(dest, exist_ok=True)
    archive = _fetch(f"{EXOMISER_DATA_URL}/{version}_phenotype.zip", annotations_dir)
    _extract_zip(archive, dest)
    if not keep:
        os.remove(archive)
    _chmod_777(os.path.join(annotations_dir, "Annotations_Exomiser"))


def install_exomiser_jar(annotations_dir: str, version: str = EXOMISER_JAR_VERSION,
                         jar_file: Optional[str] = None) -> str:
    """install-exomiser-resources: the rest-prioritiser jar in <annotationsDir>/jar/.

    ``jar_file`` (Makefile EX_RP_FILE): a jar already downloaded, symlinked
    instead of downloaded. (The Makefile symlink used an undefined variable
    for its name; the link is named after the jar here.)
    """
    jar_dir = os.path.join(annotations_dir, "jar")
    os.makedirs(jar_dir, exist_ok=True)
    if jar_file:
        _log("Custom Exomiser rest-prioritiser path provided; creating symlink")
        link = os.path.join(jar_dir, os.path.basename(jar_file))
        if os.path.lexists(link):
            os.remove(link)
        os.symlink(os.path.abspath(jar_file), link)
        return link
    name = f"exomiser-rest-prioritiser-{version}.jar"
    target = os.path.join(jar_dir, name)
    if not os.path.exists(target):
        part = target + ".part"
        download(f"{EXOMISER_JAR_URL}/{version}/{name}", part)
        os.replace(part, target)
        os.chmod(target, 0o755)
    return target


def fix_variantconvert_configs(config_dir: str) -> None:
    """Post-install edits of the variantconvert configs (install-variantconvert):
    hs1 renamed CHM13, no "chr" prefix in the CHM13 contig lines, and the empty
    ``*.local.json`` files used by the Conda recipe."""
    hs1, chm13 = os.path.join(config_dir, "hs1"), os.path.join(config_dir, "CHM13")
    if os.path.isdir(hs1) and not os.path.exists(chm13):
        os.rename(hs1, chm13)
    for path in glob.glob(os.path.join(chm13, "annotsv*.json")):
        with open(path) as f:
            text = f.read()
        new = text.replace('"##contig=<ID=chr', '"##contig=<ID=')
        if new != text:
            with open(path, "w") as f:
                f.write(new)
    for genome in ("GRCh37", "GRCh38", "CHM13"):
        for source in ("bed", "vcf"):
            for kind in ("combined", "full", "fullsplit"):
                p = os.path.join(config_dir, genome, f"annotsv3_from_{source}.{kind}.local.json")
                os.makedirs(os.path.dirname(p), exist_ok=True)
                open(p, "a").close()
    _chmod_777(config_dir)


def install_variantconvert(install_dir: str, version: str = VC_VERSION) -> bool:
    """install-variantconvert: clone into share/python3/variantconvert and pip install it."""
    _log("\nvariantconvert installation\n---------------------------")
    dest = os.path.join(install_dir, "share", "python3", "variantconvert")
    if os.path.isdir(dest):
        _log("variantconvert directory found; purging locally before re-installing.")
        shutil.rmtree(dest)
    subprocess.run(["git", "clone", "--branch", version, VARIANTCONVERT_GIT, dest], check=True)
    r = subprocess.run([sys.executable, "-m", "pip", "install", "-e", dest], capture_output=True, text=True)
    if r.returncode != 0:
        _log(r.stdout + r.stderr)
        _log("variantconvert not installed")
        return False
    flag = os.path.join(dest, "pipinstall.flag")   # checked by AnnotSV at run time
    open(flag, "a").close()
    os.chmod(flag, 0o777)
    fix_variantconvert_configs(os.path.join(dest, "src", "variantconvert", "configs"))
    _log("variantconvert installed")
    return True


TARGETS = ("human", "mouse", "exomiser", "exomiser-jar", "variantconvert")


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(prog="annotsv-install", description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("targets", nargs="+", choices=TARGETS, metavar="TARGET",
                   help="one or more of: " + ", ".join(TARGETS))
    p.add_argument("--install-dir", default=None,
                   help="AnnotSV install directory (default: $ANNOTSV, the source tree or the Python prefix)")
    p.add_argument("--annotations-dir", default=None,
                   help="where annotations go (default: <install-dir>/share/AnnotSV); "
                        "pass the same directory to annotsv --annotations-dir")
    p.add_argument("--human-version", default=HUMAN_VERSION)
    p.add_argument("--mouse-version", default=MOUSE_VERSION)
    p.add_argument("--exomiser-version", default=EXOMISER_VERSION)
    p.add_argument("--exomiser-jar-version", default=EXOMISER_JAR_VERSION)
    p.add_argument("--exomiser-jar-file", default=None,
                   help="already downloaded rest-prioritiser jar to symlink instead of downloading")
    p.add_argument("--variantconvert-version", default=VC_VERSION)
    p.add_argument("--keep-archives", action="store_true", help="keep the downloaded archives")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    a = parse_args(argv)
    install_dir = os.path.abspath(a.install_dir or default_install_dir())
    ann = os.path.abspath(a.annotations_dir or os.path.join(install_dir, "share", "AnnotSV"))
    _log(f"AnnotSV install directory: {install_dir}\nAnnotations directory: {ann}")
    targets: List[str] = list(dict.fromkeys(a.targets))
    try:
        for t in targets:
            if t == "human":
                install_organism(ann, "Human", a.human_version, a.keep_archives)
                if "exomiser" not in targets:
                    install_exomiser(ann, a.exomiser_version, a.keep_archives)
            elif t == "mouse":
                install_organism(ann, "Mouse", a.mouse_version, a.keep_archives)
            elif t == "exomiser":
                install_exomiser(ann, a.exomiser_version, a.keep_archives)
            elif t == "exomiser-jar":
                _log(f"\n--> {install_exomiser_jar(ann, a.exomiser_jar_version, a.exomiser_jar_file)}")
            elif t == "variantconvert":
                if not install_variantconvert(install_dir, a.variantconvert_version):
                    return 1
    except (OSError, urllib.error.URLError, tarfile.TarError, zipfile.BadZipFile,
            subprocess.CalledProcessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    if ann != os.path.join(install_dir, "share", "AnnotSV"):
        _log(f"\nInstallation completed. Use: annotsv --annotations-dir {ann} ...")
    else:
        _log("\nInstallation completed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
