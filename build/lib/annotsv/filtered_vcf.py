"""Port of ``share/tcl/AnnotSV/AnnotSV-filteredVCF.tcl``.

Compound heterozygosity: for each sample of ``-candidateSnvIndelSamples``,
list the heterozygous SNV/indel (``${chrom}_$pos``) located in the gene that
the SV overlaps. The Tcl globals (``L_htzPos``, ``L_samplesFromTheHeader``,
``g_AnnotSV(filteredVCFparsing)``) are the state of ``FilteredVCF``.
"""

from __future__ import annotations

import datetime
import glob
import gzip
import itertools
import re
from typing import Dict, List, Optional, Sequence

from . import tcl
from .vcf import dichotomy_search


def _now() -> str:
    return datetime.datetime.now().strftime("%B %d %Y - %H:%M")


def _as_list(value) -> List[str]:
    if isinstance(value, (list, tuple)):
        return list(value)
    return tcl.split_list(value or "")


class FilteredVCF:
    """filteredVCFannotation, with its globals as instance state."""

    def __init__(self, config) -> None:
        self.config = config
        self.parsing: Optional[str] = None     # g_AnnotSV(filteredVCFparsing)
        self.htz_pos: Dict[object, List[str]] = {}
        self.samples_from_header: Optional[List[str]] = None

    def _files(self) -> List[str]:
        # Tcl: [eval glob -nocomplain $g_AnnotSV(candidateSnvIndelFiles)]
        files = []
        for item in _as_list(self.config["candidateSnvIndelFiles"]):
            matches = glob.glob(item)
            files.extend(matches if matches else [])
        return files

    def _parse(self) -> None:
        c = self.config
        self.parsing = "done"
        samples = _as_list(c["candidateSnvIndelSamples"])
        print(f'...parsing of candidateSnvIndelFiles for "{" ".join(samples)}" ({_now()})')
        self.htz_pos["FULL"] = [""] * len(samples)
        gt_absent = True
        sv_min_size = tcl.to_number(str(c["SVminSize"]))
        pass_only = str(c["snvIndelPASS"]) == "1"
        samples_ok: List[str] = []
        i_sample: Dict[str, int] = {}
        for vcf in self._files():
            print(f"\t...parse all positions from {vcf} ({_now()})")
            opener = gzip.open if re.search(r".gz$", vcf, re.I) else open
            with opener(vcf, "rt") as f:
                for raw in f:
                    line = raw[:-1] if raw.endswith("\n") else raw
                    if line[:6] == "#CHROM":
                        cols = line.split("\t")
                        samples_ok = []
                        for sample in samples:
                            i_sample[sample] = cols.index(sample) if sample in cols else -1
                            if i_sample[sample] != -1:
                                samples_ok.append(sample)
                        continue
                    if line.startswith("#") or line == "":
                        continue
                    ls = line.split("\t")
                    get = lambda j: ls[j] if j < len(ls) else ""  # noqa: E731
                    chrom = re.sub("chr", "", get(0), count=1, flags=re.I)
                    pos = get(1)
                    filt = get(6)
                    fmt = get(8)
                    ref, alt = get(3), get(4)
                    # Only SNV/indel (skip SVs)
                    if re.search(r"<|\[|\]", alt):
                        continue
                    if abs(len(ref) - len(alt)) > sv_min_size:
                        continue
                    if pass_only and filt != "PASS":
                        continue
                    for sample in samples_ok:
                        sample_data = get(i_sample[sample])
                        fmt_fields = fmt.split(":") if fmt else []
                        if "GT" not in fmt_fields:
                            continue
                        j_gt = fmt_fields.index("GT")
                        gt_absent = False
                        parts = sample_data.split(":") if sample_data else []
                        gt = parts[j_gt] if j_gt < len(parts) else ""
                        alleles = tcl.split_chars(gt, "[/|]") if gt else []
                        a0 = alleles[0] if len(alleles) > 0 else ""
                        a1 = alleles[1] if len(alleles) > 1 else ""
                        if a0 != a1 and gt != "":
                            self.htz_pos.setdefault((chrom, sample), []).append(pos)
        if gt_absent:
            print("\t...WARNING: The SNV/indel genotype is not indicated in the FORMAT column under the “GT” field")
            print("\t            => Compound heterozygosity analysis won't be processed!")
            c["candidateSnvIndelFiles"] = []
            c["candidateSnvIndelSamples"] = []
        print("")

    def annotate(self, gene_chrom: str, gene_start="", gene_end="", line: Sequence[str] = (),
                 header_output: str = "") -> Optional[str]:
        """filteredVCFannotation. ``line`` is the output line as a list.

        Returns the tab-joined compound-htz columns, or ``None`` when the SV
        input is a BED (the analysis is then switched off).
        """
        c = self.config
        if re.search(r".bed(.gz)?", c["SVinputFile"], re.I):
            print("\t...WARNING: The compound heterozygosity analysis is not available from a BED SV input file")
            print("\t            => That analysis won't be processed!\n")
            c["candidateSnvIndelFiles"] = []
            c["candidateSnvIndelSamples"] = []
            self.parsing = "useless"
            return None

        if self.parsing is None:
            self._parse()

        if gene_chrom == "FULL":
            return "\t".join(self.htz_pos["FULL"])

        if self.samples_from_header is None:
            samples = re.sub(r".*FORMAT\t", "", header_output, count=1, flags=re.S)
            samples = re.sub(r"\tAnnotation_mode.*", "", samples, count=1, flags=re.S)
            self.samples_from_header = samples.split("\t")

        detected = set()
        fcol = c.get("formatColNumber")
        if fcol is not None:
            fcol = int(fcol)
            fmt = line[fcol] if fcol < len(line) else ""
            fmt_fields = fmt.split(":") if fmt else []
            i_gt = fmt_fields.index("GT") if "GT" in fmt_fields else -1
            values = list(line[fcol + 1:len(line) - 13]) if len(line) - 13 > fcol + 1 else []
            for sample, value in itertools.zip_longest(self.samples_from_header, values, fillvalue=""):
                parts = value.split(":") if value else []
                gt = parts[i_gt] if 0 <= i_gt < len(parts) else ""
                if re.search(r"[1-9]", gt):
                    detected.add(sample)

        out: List[str] = []
        for sample in _as_list(c["candidateSnvIndelSamples"]):
            if sample not in detected or (gene_chrom, sample) not in self.htz_pos:
                out.append("")
                continue
            positions = self.htz_pos[(gene_chrom, sample)]
            i_first = dichotomy_search(gene_start, positions)
            i_last = dichotomy_search(gene_end, positions)
            if str(gene_start) != (positions[i_first] if 0 <= i_first < len(positions) else ""):
                i_first += 1
            lo = max(i_first, 0)
            chunk = positions[lo:i_last + 1] if i_last >= lo else []
            out.append(";".join(f"{gene_chrom}_{p}" for p in sorted(set(chunk))))
        return "\t".join(out)
