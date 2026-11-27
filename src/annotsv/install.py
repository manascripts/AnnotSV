"""Port of the annotation installer (``bin/INSTALL_annotations.sh``) and of the
data installation targets of the ``Makefile``.

The code itself is installed with ``pip install .`` (it replaces the
``install-executable``, ``install-tcl-toolbox``, ``install-configfile``... targets).
This command installs what pip cannot: the annotation bundles and the external
tools AnnotSV downloads.

    annotsv-install [--annotations-dir DIR] [TARGET ...]

Targets (Makefile target in brackets):

* (none)/``base``     variantconvert + the Exomiser jar     [install]
* ``human``           Exomiser phenotype data + Annotations_Human_<v>.tar.gz
                      [install-human-annotation, which also runs install-exomiser]
* ``mouse``           Annotations_Mouse_<v>.tar.gz          [install-mouse-annotation]
* ``all-annotations`` human + mouse                         [install-all-annotations]
* ``exomiser``        <v>_phenotype.zip only                [install-exomiser]
* ``exomiser-jar``    exomiser-rest-prioritiser jar          [install-exomiser-resources]
* ``variantconvert``  git clone + pip install of variantconvert [install-variantconvert]

As with ``make``, the annotations are never installed implicitly:
``annotsv-install base human`` is ``make install install-human-annotation``.
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
import time
import urllib.error
import urllib.request
import zipfile
from typing import IO, List, Optional, Sequence

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


# Progress display: "auto" (a redrawn bar on a terminal, a line every 10% in a
# log file), "bar", "log" or "none". Set by --progress.
PROGRESS = "auto"


def _size(n: float) -> str:
    for unit in ("B", "KiB", "MiB", "GiB"):
        if n < 1024 or unit == "GiB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024
    return ""


def _duration(sec: float) -> str:
    sec = int(sec)
    h, rem = divmod(sec, 3600)
    m, s = divmod(rem, 60)
    return f"{h}h{m:02d}m" if h else f"{m}m{s:02d}s" if m else f"{s}s"


class Progress:
    """Download progress, the stdlib equivalent of curl's progress meter.

    On a terminal the bar is redrawn in place (at most 10 times a second).
    Otherwise (cluster job logs, ``| tee``...) one line is written per 10%,
    or once a minute when the size is unknown, so logs stay short.
    """

    WIDTH = 30

    def __init__(self, name: str, total: Optional[int], start: int = 0,
                 mode: Optional[str] = None, stream: Optional[IO[str]] = None):
        self.name, self.total, self.done, self.start = name, total, start, start
        self.stream = stream or sys.stdout
        mode = mode or PROGRESS
        if mode == "auto":
            isatty = getattr(self.stream, "isatty", None)
            mode = "bar" if isatty and isatty() else "log"
        self.mode = mode
        self.t0 = self._last = time.monotonic()
        self._next_pct = 10
        if start and total:
            self._next_pct = (start * 100 // total) // 10 * 10 + 10
        if self.mode == "log":
            what = f"{_size(total)}" if total else "unknown size"
            resumed = f", resuming at {_size(start)}" if start else ""
            self._write(f"Downloading {name} ({what}{resumed})\n")

    def _write(self, text: str) -> None:
        self.stream.write(text)
        self.stream.flush()

    def _rate(self, now: float) -> float:
        return (self.done - self.start) / max(now - self.t0, 1e-6)

    def _line(self, now: float) -> str:
        rate = self._rate(now)
        if self.total:
            frac = min(self.done / self.total, 1.0)
            eta = (self.total - self.done) / rate if rate > 0 else 0
            return (f"{frac * 100:3.0f}%  {_size(self.done)} / {_size(self.total)}"
                    f"  {_size(rate)}/s  ETA {_duration(eta)}")
        return f"{_size(self.done)}  {_size(rate)}/s"

    def update(self, n: int) -> None:
        self.done += n
        if self.mode == "none":
            return
        now = time.monotonic()
        if self.mode == "bar":
            if now - self._last >= 0.1:
                self._last = now
                self._draw(now)
        elif self.total:
            pct = self.done * 100 // self.total
            if pct >= self._next_pct and pct < 100:
                self._next_pct = pct // 10 * 10 + 10
                self._write(f"  {self.name}: {self._line(now)}\n")
        elif now - self._last >= 60:
            self._last = now
            self._write(f"  {self.name}: {self._line(now)}\n")

    def _draw(self, now: float) -> None:
        if self.total:
            filled = int(self.WIDTH * min(self.done / self.total, 1.0))
            bar = "[" + "#" * filled + "-" * (self.WIDTH - filled) + "] "
        else:
            bar = ""
        self._write(f"\r{self.name}  {bar}{self._line(now)}\033[K")

    def close(self) -> None:
        if self.mode == "none":
            return
        now = time.monotonic()
        if self.mode == "bar":
            self._draw(now)
            self._write("\n")
        else:
            took = _duration(now - self.t0)
            self._write(f"  {self.name}: done, {_size(self.done)} in {took}"
                        f" ({_size(self._rate(now))}/s)\n")


def download(url: str, dest: str, progress: Optional[str] = None) -> str:
    """``curl -C - -LO url``: download into ``dest``, resuming a partial file,
    with a progress display (``progress``: see ``PROGRESS``)."""
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
        resumed = bool(done) and getattr(resp, "status", 200) == 206
        start = done if resumed else 0
        length = resp.headers.get("Content-Length")
        total = start + int(length) if length and length.isdigit() else None
        bar = Progress(url.rsplit("/", 1)[-1], total, start, progress)
        with open(dest, "ab" if resumed else "wb") as out:
            while True:
                chunk = resp.read(1 << 20)
                if not chunk:
                    break
                out.write(chunk)
                bar.update(len(chunk))
        bar.close()
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
    _log(f"Extracting {os.path.basename(archive)}")
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
    _log(f"Extracting {os.path.basename(archive)}")
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
    _log("\nInstallation of the Exomiser rest-prioritiser:\n")
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
    if os.path.exists(target):
        _log(f"{name} already installed")
    else:
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


TARGETS = ("base", "human", "mouse", "all-annotations", "exomiser", "exomiser-jar", "variantconvert")

# "make install" installs variantconvert then the Exomiser resources; the
# annotations are separate targets that "make install-all-annotations" groups.
ALIASES = {
    "base": ["variantconvert", "exomiser-jar"],
    "all-annotations": ["human", "mouse"],
}


def expand_targets(targets: Sequence[str]) -> List[str]:
    """Targets in the order make runs them, each once. No target is ``make install``."""
    out: List[str] = []
    for t in targets or ["base"]:
        for x in ALIASES.get(t, [t]):
            if x == "human":
                # install-human-annotation depends on install-exomiser: run first
                out.append("exomiser")
            out.append(x)
    return list(dict.fromkeys(out))


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(prog="annotsv-install", description=__doc__.split("\n\n")[0],
                                epilog=__doc__.split("\n\n", 3)[3],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    # choices= is checked by hand: with nargs="*", Python < 3.12 rejects no target
    p.add_argument("targets", nargs="*", metavar="TARGET",
                   help="zero or more of: " + ", ".join(TARGETS) + " (default: base)")
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
    p.add_argument("--progress", choices=("auto", "bar", "log", "none"), default="auto",
                   help="download progress: a bar on a terminal, a line every 10%% in a log "
                        "file (auto, the default), or none")
    a = p.parse_args(argv)
    for t in a.targets:
        if t not in TARGETS:
            p.error(f"argument TARGET: invalid choice: {t!r} (choose from {', '.join(TARGETS)})")
    return a


def main(argv: Optional[Sequence[str]] = None) -> int:
    global PROGRESS
    a = parse_args(argv)
    PROGRESS = a.progress
    install_dir = os.path.abspath(a.install_dir or default_install_dir())
    ann = os.path.abspath(a.annotations_dir or os.path.join(install_dir, "share", "AnnotSV"))
    _log(f"AnnotSV install directory: {install_dir}\nAnnotations directory: {ann}")
    targets = expand_targets(a.targets)
    annotations = {"human", "mouse"} & set(targets)
    vc_failed = False
    try:
        for t in targets:
            if t == "human":
                install_organism(ann, "Human", a.human_version, a.keep_archives)
            elif t == "mouse":
                install_organism(ann, "Mouse", a.mouse_version, a.keep_archives)
            elif t == "exomiser":
                install_exomiser(ann, a.exomiser_version, a.keep_archives)
            elif t == "exomiser-jar":
                _log(f"--> {install_exomiser_jar(ann, a.exomiser_jar_version, a.exomiser_jar_file)}")
            elif t == "variantconvert":
                # As with make, a failed pip install does not stop the installation
                vc_failed = not install_variantconvert(install_dir, a.variantconvert_version)
    except (OSError, urllib.error.URLError, tarfile.TarError, zipfile.BadZipFile,
            subprocess.CalledProcessError) as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        return 1
    _log("\nDone")
    if vc_failed:
        _log("\nWARNING: variantconvert not installed (no VCF output); see the pip output above.")
    custom = ann != os.path.join(install_dir, "share", "AnnotSV")
    opt = f" --annotations-dir {ann}" if custom else ""
    if not annotations:
        # install-done
        _log("\nWARNING: Annotations need to be installed:")
        _log(f"annotsv-install{opt} human")
        _log(f"annotsv-install{opt} mouse")
    elif custom:
        _log(f"\nInstallation completed. Use: annotsv --annotations-dir {ann} ...")
    else:
        _log("\nInstallation completed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
