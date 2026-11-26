"""Reading of the SV input file (VCF or BED).

Port of the input step of ``bin/AnnotSV`` and of the ``createBEDinputHeaderFile``
and ``addNAinSamplesIDbedCol`` procs from ``AnnotSV-config.tcl``.

The Tcl code kept its results in globals; here they go to the config
(``bedFile``, ``NAbedFile``, ``samplesidBEDcol``) and to the returned
``SVInput`` (``VCFheader``, ``g_SVLEN``, ``headerFileToRemove``).
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from typing import Dict, Optional

from . import tcl
from .general import first_line_from_file, write_text_in_file
from .vcf import vcfs_to_bed

NA_CHUNK = 100000


@dataclass
class SVInput:
    bed_file: str
    # VCF input only (the Tcl "VCFheader" and "g_SVLEN" globals)
    vcf_header: Optional[str] = None
    svlen: Dict[str, object] = field(default_factory=dict)
    # BED input only
    header_file: Optional[str] = None       # the *.header.tsv written in outputDir
    header_file_to_remove: bool = False
    add_na_in_bed: bool = False
    unannotated_file: str = ""


def _header_tsv_name(bed_file: str, suffix: str = ".header.tsv") -> str:
    # Tcl: regsub -nocase ".bed$" ... (the dot matches any character)
    return re.sub(r".bed$", lambda _m: suffix, bed_file, count=1, flags=re.I)


def bed_input_header_path(config) -> str:
    return f"{config['outputDir']}/{os.path.basename(_header_tsv_name(config['bedFile']))}"


def create_bed_input_header_file(config) -> SVInput:
    """createBEDinputHeaderFile: write ``<outputDir>/<input>.header.tsv``.

    Uses the last "#" line before the first data line as the header, or builds
    ``SV_chrom SV_start SV_end [SV_type|Samples_ID|user#k ...]`` when there is
    none. Adds a Samples_ID column when no column holds the sample IDs.
    """
    c = config
    result = SVInput(bed_file=c["bedFile"])
    header_file = bed_input_header_path(c)
    header: Optional[str] = None
    with open(c["bedFile"], "r") as f:
        for raw in f:
            line = raw[:-1] if raw.endswith("\n") else raw
            if line == "":
                continue
            if line.startswith("#"):
                header = line
                continue
            break

    if header is not None:
        if tcl.compare(c["samplesidBEDcol"], "-1") == 0:
            # Tcl: [lsearch -exact "Samples_ID" [split $header "\t"]] - the
            # arguments are swapped, so it only "finds" a header that is
            # exactly "Samples_ID" (index 0). Kept as is.
            i_samplesid = 0 if tcl.format_list(header.split("\t")) == "Samples_ID" else -1
            if i_samplesid != -1:
                c["samplesidBEDcol"] = str(i_samplesid + 1)
            else:
                header += "\tSamples_ID"
                result.add_na_in_bed = True
                c["samplesidBEDcol"] = str(len(header.split("\t")))
    else:
        header = "SV_chrom\tSV_start\tSV_end"
        bed_length = len(first_line_from_file(c["bedFile"]).split("\t"))
        k = 1
        for i in range(4, bed_length + 1):
            if str(i) == str(c["svtBEDcol"]):
                header += "\tSV_type"
            elif str(i) == str(c["samplesidBEDcol"]):
                header += "\tSamples_ID"
            else:
                header += f"\tuser#{k}"
                k += 1
        if tcl.compare(c["samplesidBEDcol"], "-1") == 0:
            header += "\tSamples_ID"
            result.add_na_in_bed = True
            c["samplesidBEDcol"] = str(bed_length + 1)
    write_text_in_file(header, header_file)
    result.header_file = header_file
    result.header_file_to_remove = True
    return result


def add_na_in_samples_id_bed_col(config, state: SVInput) -> None:
    """addNAinSamplesIDbedCol: when the BED has no Samples_ID column, write
    ``<input>.NA.bed`` (next to the input) with an extra "NA" column, rename the
    header file to ``*.NA.header.tsv`` and use the new BED from now on."""
    if not state.add_na_in_bed:
        return
    c = config
    na_bed = _header_tsv_name(c["bedFile"], ".NA.bed")
    c["NAbedFile"] = na_bed
    if os.path.exists(na_bed):
        os.remove(na_bed)
    to_write = []
    i = 0
    with open(c["bedFile"], "r") as f:
        for raw in f:
            line = raw[:-1] if raw.endswith("\n") else raw
            if line == "" or line.startswith("#"):
                continue
            to_write.append(f"{line}\tNA")
            i += 1
            if i > NA_CHUNK:
                write_text_in_file("\n".join(to_write), na_bed)
                to_write = []
                i = 0
    write_text_in_file("\n".join(to_write), na_bed)

    header_file = bed_input_header_path(c)
    na_header_file = f"{c['outputDir']}/{os.path.basename(_header_tsv_name(c['bedFile'], '.NA.header.tsv'))}"
    os.replace(header_file, na_header_file)
    state.header_file = na_header_file
    c["bedFile"] = na_bed
    state.bed_file = na_bed


def read_sv_input(config, now=None) -> SVInput:
    """The "VCF or BED input" step of ``bin/AnnotSV``.

    VCF: converted to a BED with ``vcfs_to_bed``. BED: a header file is written
    in outputDir and a Samples_ID column of "NA" is added when missing.
    """
    c = config
    if re.search(r"\.vcf(.gz)?$", c["SVinputFile"], re.I):
        res = vcfs_to_bed(c, c["SVinputFile"], now=now)
        c["bedFile"] = res.bed_file
        return SVInput(bed_file=res.bed_file, vcf_header=res.vcf_header, svlen=res.svlen,
                       unannotated_file=res.unannotated_file)

    c["bedFile"] = c["SVinputFile"]
    stale_header = _header_tsv_name(c["bedFile"])
    if c["bedFile"] != stale_header and os.path.exists(stale_header):
        os.remove(stale_header)
    state = create_bed_input_header_file(c)
    add_na_in_samples_id_bed_col(c, state)
    return state
