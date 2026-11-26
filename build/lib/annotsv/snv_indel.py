"""Port of ``VCFannotation`` from ``share/tcl/AnnotSV/AnnotSV-vcf.tcl``.

For each sample of ``-snvIndelSamples``, count the homozygous and
heterozygous SNV/indel of ``-snvIndelFiles`` located in an SV, and the total
count in the cohort. The values are reported for deletions only.

The Tcl globals (``lPos``, ``L_allSamples``, ``g_AnnotSV(vcfParsing)``) are
the state of ``SnvIndelCounts``. Tcl quirks kept on purpose:

* the positions are kept in the VCF order (``DichotomySearch`` expects them
  sorted, as they are for a sorted VCF);
* a ``./.`` genotype counts as homozygous, ``0/.`` as heterozygous;
* the SV coordinates are compared to the VCF positions as they are (BED
  start against 1-based POS).
"""

from __future__ import annotations

import datetime
import glob
import gzip
import os
import re
import shutil
import subprocess
import time
from typing import Dict, List, Optional

from . import tcl
from .general import normalize_svtype, write_text_in_file
from .vcf import dichotomy_search


def _now() -> str:
    return datetime.datetime.now().strftime("%B %d %Y - %H:%M")


def _as_list(value) -> List[str]:
    if isinstance(value, (list, tuple)):
        return list(value)
    return tcl.split_list(value or "")


def _gzip_file(path: str) -> None:
    """``exec gzip $path``: replaces ``path`` by ``path.gz``."""
    with open(path, "rb") as src, gzip.open(path + ".gz", "wb") as dst:
        shutil.copyfileobj(src, dst)
    os.remove(path)


def _gzip_copy(src_path: str, dst_path: str, mode: str = "wb") -> None:
    """``cat $src | gzip > $dst`` (``>>`` with mode "ab")."""
    with open(src_path, "rb") as src, gzip.open(dst_path, mode) as dst:
        shutil.copyfileobj(src, dst)


def _gz_lines(path: str):
    with gzip.open(path, "rt", newline="") as f:
        for raw in f:
            yield raw[:-1] if raw.endswith("\n") else raw


def _norm_message(c, src: str, dst: str) -> str:
    """``catch {eval exec $bcftools norm -m -both $src > $dst} Message``: the message."""
    cmd = tcl.split_list(str(c.get("bcftools", "bcftools"))) + ["norm", "-m", "-both", src]
    with open(dst, "wb") as out:
        try:
            r = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE)
        except OSError:
            return f"couldn't execute \"{cmd[0]}\": no such file or directory"
    err = r.stderr.decode(errors="replace")
    err = err[:-1] if err.endswith("\n") else err
    if err:
        return err
    return "child process exited abnormally" if r.returncode else ""


class SnvIndelCounts:
    """VCFannotation, with its globals as instance state."""

    def __init__(self, config) -> None:
        self.config = config
        self.parsed = False                         # g_AnnotSV(vcfParsing)
        self.pos: Dict[tuple, List[str]] = {}       # lPos($chrom,$GT,$sample)
        self.all_samples: List[str] = []            # L_allSamples

    # ------------------------------------------------------------------
    def _files(self) -> List[str]:
        files: List[str] = []
        for item in _as_list(self.config["snvIndelFiles"]):
            files.extend(glob.glob(item))
        return files

    def _parse(self) -> None:
        c = self.config
        self.parsed = True
        out_dir = c["outputDir"]
        print(f'\n\n...parsing of snvIndelFiles for "{" ".join(_as_list(c["snvIndelSamples"]))}" ({_now()})')
        sv_min_size = tcl.to_number(str(c["SVminSize"]))
        pass_only = str(c.get("snvIndelPASS", "0")) == "1"
        bedtools = tcl.split_list(str(c.get("bedtools", "bedtools")))
        bcftools = str(c.get("bcftools", "bcftools"))

        for vcf in self._files():
            # Split multiallelic sites into multiple rows
            print(f"\t...split multiallelic sites into multiple rows for {vcf} ({_now()})")
            stamp = int(time.time())
            tmp1 = f"{out_dir}/{os.path.basename(vcf)}.{stamp}.separate.vcf"
            tmp_gz1 = f"{out_dir}/{os.path.basename(vcf)}.{stamp}.separate.vcf.gz"
            message = _norm_message(c, vcf, tmp1)
            failed = re.search("error|fail", message, re.I) is not None
            if os.path.getsize(tmp1) == 0 or failed:
                print("\t   -- VCFannotation --")
                print(f"\t   {bcftools} norm -m -both {vcf} > {tmp1}")
                if os.path.getsize(tmp1) == 0:
                    print(f"\t   {tmp1}: file empty.")
                for line in message.split("\n"):
                    if re.search("error|fail", line, re.I):
                        print(f"\t   {line}")
                print("\t   => No multiallelic treatment done.")
                if re.search(".gz$", vcf, re.I):
                    shutil.copyfile(vcf, tmp_gz1)
                else:
                    _gzip_copy(vcf, tmp_gz1)
            else:
                if os.path.exists(tmp1 + ".gz"):
                    os.remove(tmp1 + ".gz")
                _gzip_file(tmp1)
            if os.path.exists(tmp1):
                os.remove(tmp1)

            # Parsing of $vcfF
            print(f"\t...parsing of {vcf} ({_now()})")
            i_intersect = i_sv = i_loaded = i_not_pass = i_not_gt = 0

            # VCF should not contain "chr" prefix
            tmp_gz2 = re.sub(".vcf.gz$", lambda _m: ".withoutChr.vcf.gz", tmp_gz1, count=1)
            if any(line.startswith("chr") for line in _gz_lines(tmp_gz1)):
                with gzip.open(tmp_gz2, "wb") as dst:
                    for raw in gzip.open(tmp_gz1, "rb"):
                        dst.write(raw[3:] if raw.startswith(b"chr") else raw)
            else:
                tmp_gz2 = tmp_gz1

            # Header for the intersection with the SV bedfile
            tmp_gz3 = re.sub(".vcf.gz$", lambda _m: ".intersect.vcf.gz", tmp_gz2, count=1)
            tmp3 = re.sub(".gz$", "", tmp_gz3, count=1)
            header = ""
            for line in _gz_lines(tmp_gz2):
                if line.startswith("#"):
                    header = line
                else:
                    break
            self.all_samples.extend(tcl.split_list(header)[9:])
            write_text_in_file(header, tmp3)
            if os.path.exists(tmp_gz3):
                os.remove(tmp_gz3)
            _gzip_file(tmp3)

            # Intersect with the bedfile
            with open(tmp3, "wb") as out:
                r = subprocess.Popen(bedtools + ["intersect", "-a", "stdin", "-b", c["bedFile"]],
                                     stdin=subprocess.PIPE, stdout=out, stderr=subprocess.PIPE)
                with gzip.open(tmp_gz2, "rb") as src:
                    _, err = r.communicate(src.read())
            err_text = err.decode(errors="replace").rstrip("\n")
            if re.search("error", err_text, re.I):
                write_text_in_file(err_text, tmp3 + ".error")
                print("-- VCFannotation --")
                print(f"gunzip -c {tmp_gz2} | {' '.join(bedtools)} intersect -a stdin -b {c['bedFile']} > {tmp3}")
                print(f"Error: look at {tmp3}.error")

            for path in (tmp_gz1, tmp_gz2):
                if os.path.exists(path):
                    os.remove(path)

            if os.path.getsize(tmp3) == 0:
                print("\t\t-> no variant in the intersection")
                os.remove(tmp3)
                os.remove(tmp_gz3)
                continue
            _gzip_copy(tmp3, tmp_gz3, "ab")
            os.remove(tmp3)

            # load SNV/Indel in memory
            samples_1vcf: List[str] = []
            i_sample: Dict[str, int] = {}
            for line in _gz_lines(tmp_gz3):
                ls = line.split("\t")
                if line[:6] == "#CHROM":
                    samples_1vcf = ls[9:]
                    for k, sample in enumerate(samples_1vcf, start=9):
                        i_sample[sample] = k
                    continue
                if line.startswith("#") or line == "":
                    continue
                i_intersect += 1
                get = lambda j: ls[j] if 0 <= j < len(ls) else ""  # noqa: E731
                chrom = re.sub("chr", "", get(0), count=1, flags=re.I)
                pos = get(1)
                ref, alt = get(3), get(4)
                if re.search(r"<|\[|\]", alt):
                    i_sv += 1
                    continue
                if abs(len(ref) - len(alt)) > sv_min_size:
                    i_sv += 1
                    continue
                filt = get(6)
                fmt = get(8).split(":")
                if "GT" not in fmt:
                    i_not_gt += 1
                    continue
                j_gt = fmt.index("GT")
                if pass_only and filt != "PASS":
                    i_not_pass += 1
                    continue
                i_loaded += 1
                for sample in samples_1vcf:
                    data = get(i_sample[sample]).split(":")
                    gt = tcl.split_chars(data[j_gt] if j_gt < len(data) else "", "/|\\")
                    a = gt[0] if len(gt) > 0 else ""
                    b = gt[1] if len(gt) > 1 else ""
                    if a == "" or b == "":
                        continue
                    if a == "0" and b == "0":
                        continue
                    kind = "htz" if a != b else "hom"
                    self.pos.setdefault((chrom, kind, sample), []).append(pos)
            os.remove(tmp_gz3)
            if i_intersect:
                print(f"\t\t-> {i_intersect} SNV/indel located in all the SV (potentially with redundancy "
                      "depending on SV overlap)")
            if i_loaded:
                print(f"\t\t-> {i_loaded} SNV/indel loaded")
            if i_sv:
                print(f"\t\t-> {i_sv} SV excluded (considering only SNV/indel from the VCF)")
            if i_not_pass:
                print(f"\t\t-> {i_not_pass} SNV/indel excluded because of the FILTER value not equal to \"PASS\"")
            if i_not_gt:
                print(f"\t\t-> {i_not_gt} SNV/indel excluded because of the absence of GT information")

    # ------------------------------------------------------------------
    def _in_sv(self, key, start, end) -> List[str]:
        positions = self.pos[key]
        i_first = dichotomy_search(start, positions)
        i_last = dichotomy_search(end, positions)
        if str(start) != (positions[i_first] if 0 <= i_first < len(positions) else ""):
            i_first += 1
        # lrange: a negative first index is 0
        lo = max(i_first, 0)
        return sorted(set(positions[lo:i_last + 1])) if i_last >= lo else []

    def header(self) -> str:
        text = ""
        for sample in _as_list(self.config["snvIndelSamples"]):
            text += (f"\tCount_hom({sample})\tCount_htz({sample})\tCount_htz/allHom({sample})"
                     "\tCount_htz/total(cohort)")
        return text + "\tCount_total(cohort)"

    def text(self, chrom: str, start, end, svtype: str) -> str:
        """VCFannotation: the tab-joined count columns."""
        if not self.parsed:
            self._parse()
        count_hom: Dict[str, int] = {}
        count_htz: Dict[str, int] = {}
        cohort = set()
        for sample in self.all_samples:
            count_hom[sample] = 0
            if (chrom, "hom", sample) in self.pos:
                found = self._in_sv((chrom, "hom", sample), start, end)
                count_hom[sample] = len(found)
                cohort.update(found)
            count_htz[sample] = 0
            if (chrom, "htz", sample) in self.pos:
                found = self._in_sv((chrom, "htz", sample), start, end)
                count_htz[sample] = len(found)
                cohort.update(found)
        total = len(cohort)

        samples = _as_list(self.config["snvIndelSamples"])
        text = ""
        if "DEL" in normalize_svtype(svtype):
            for sample in samples:
                hom, htz = count_hom[sample], count_htz[sample]  # Tcl: "no such variable" otherwise
                if total > 50:
                    all_hom = total - htz
                    htz_hom = "NA" if all_hom == 0 else "%.4f" % (htz * 1.0 / all_hom)
                    htz_tot = "%.4f" % (htz * 1.0 / total)
                    text += f"{hom}\t{htz}\t{htz_hom}\t{htz_tot}\t"
                else:
                    text += f"{hom}\t{htz}\tNA\tNA\t"
            return text + str(total)
        return "NA\tNA\tNA\tNA\t" * len(samples) + "NA"
