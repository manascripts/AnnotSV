"""Port of the SV VCF input reading from ``share/tcl/AnnotSV/AnnotSV-vcf.tcl``.

``vcfs_to_bed`` replaces the ``VCFsToBED`` proc: it turns the SV VCF input
file into the BED file that the annotation steps read, writes the variants it
skips to the ``*.unannotated.tsv`` file, and returns what the Tcl code kept in
the ``VCFheader`` and ``g_SVLEN`` globals.

``dichotomy_search`` replaces ``DichotomySearch``. ``VCFannotation`` (the
SNV/indel counts inside an SV) is ported with the SNV/indel annotation module.
"""

from __future__ import annotations

import datetime
import gzip
import os
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

from . import tcl
from .config import ConfigExit
from .general import AnnotSVIDRegistry, is_an_empty_file, normalize_svtype, write_text_in_file

ORDERED_CHROMS = "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M MT".split()
CLASSIC_HEADER = "#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT".split()
WRITE_CHUNK = 500000

_BASES = r"[acgtnACGTN]"
_SQUARE_BRACKET_ALT = re.compile(r"([acgtnACGTN]*)(\[|\])([^:]+):([0-9]+)(\[|\])([acgtnACGTN]*)")


def dichotomy_search(number, ordered: Sequence) -> int:
    """DichotomySearch: index after which ``number`` fits in ``ordered``.

    Returns -1 when ``number`` is below the first element.
    """
    if not ordered:
        # Tcl compares against "" and returns -1 (j = -1 on the second test).
        return -1
    i, j = 0, len(ordered) - 1
    if tcl.compare(str(number), str(ordered[0])) < 0:
        return -1
    if tcl.compare(str(number), str(ordered[j])) >= 0:
        return j
    k = (i + j) // 2
    if tcl.compare(str(number), str(ordered[k])) < 0:
        j = k
    else:
        i = k
    while j - i > 1:
        k = (i + j) // 2
        if tcl.compare(str(number), str(ordered[k])) < 0:
            j = k
        else:
            i = k
    return i


# --------------------------------------------------------------------------
# Small Tcl helpers
# --------------------------------------------------------------------------
def _int(value) -> int:
    """``expr`` integer arithmetic on a value read from the VCF."""
    if isinstance(value, int):
        return value
    n = tcl.to_number(value)
    if not isinstance(n, int):
        raise ValueError(f'expected integer but got "{value}"')
    return n


def _num(value):
    if isinstance(value, (int, float)):
        return value
    n = tcl.to_number(value)
    if n is None:
        raise ValueError(f'expected number but got "{value}"')
    return n


def _lt(a, b) -> bool:
    return tcl.compare(str(a), str(b)) < 0


def _idx(items: List[str], i: int) -> str:
    """``lindex`` (empty string when out of range)."""
    return items[i] if 0 <= i < len(items) else ""


def _first_sub(pattern: str, repl: str, text: str, flags: int = 0) -> str:
    """``regsub`` without ``-all`` (the replacement is literal)."""
    return re.sub(pattern, lambda _m: repl, text, count=1, flags=flags)


def _regsub_on_list(pattern: str, repl: str, items: List[str]) -> List[str]:
    """``regsub $pattern $list $repl list``: first match across the elements."""
    for n, item in enumerate(items):
        if re.search(pattern, item):
            items = list(items)
            items[n] = _first_sub(pattern, repl, item)
            return items
    return items


def _tcl_regsub_safe(pattern: str, repl: str, text: str, flags: int = 0) -> str:
    """``regsub`` with a pattern built from VCF data (REF used as a regexp)."""
    try:
        return _first_sub(pattern, repl, text, flags)
    except re.error:
        return text


def _ci_incr(value, delta):
    """``catch {incr value delta}``: unchanged when Tcl's incr would fail."""
    a, d = tcl.to_number(str(value)), tcl.to_number(str(delta))
    if isinstance(a, int) and isinstance(d, int):
        return a + d
    return value


@dataclass
class VCFtoBEDResult:
    bed_file: str
    vcf_header: str
    # g_SVLEN: AnnotSV_ID -> SV length (int, or "" / text kept as read)
    svlen: Dict[str, object] = field(default_factory=dict)
    unannotated_file: str = ""


def unannotated_output_file(config) -> str:
    path = f"{config['outputDir']}/{config['outputFile']}"
    if re.search("annotated", path):
        return _first_sub("annotated", "unannotated", path)
    return _first_sub(r".tsv$", ".unannotated.tsv", path)


def _open_vcf(path: str):
    if re.search(r".gz$", path, re.I):
        return gzip.open(path, "rt")
    return open(path, "r")


def _read_lines(f):
    """``gets`` until EOF: lines without "\\n", then the final empty read."""
    for raw in f:
        yield raw[:-1] if raw.endswith("\n") else raw
    yield ""


def vcfs_to_bed(config, vcf_file: Optional[str] = None,
                now: Optional[datetime.datetime] = None) -> VCFtoBEDResult:
    """VCFsToBED: convert the SV VCF input into ``*_AnnotSV_inputSVfile.bed``.

    BED is 0-based, half open; VCF is 1-based, inclusive. Raises ``ConfigExit``
    where the Tcl code exits (bad header: 2, no SV left: 0).
    """
    c = config
    vcf_file = vcf_file or c["SVinputFile"]
    stamp = (now or datetime.datetime.now()).strftime("%Y%m%d-%H%M%S")
    sv_bed = f"{c['outputDir']}/{stamp}_{os.getpid()}_AnnotSV_inputSVfile.bed"
    if os.path.exists(sv_bed):
        os.remove(sv_bed)
    unannotated = unannotated_output_file(c)
    if os.path.exists(unannotated):
        os.remove(unannotated)

    def skip(text: str) -> None:
        write_text_in_file(text, unannotated)

    sv_min_size = _num(c["SVminSize"])
    max_size = _num(c["bracketedAltMaxSize"]) * 1000000
    brack_mode = c["bracketedAltMode"]
    sv_input_info = str(c["SVinputInfo"]) == "1"
    include_ci = str(c["includeCI"]) == "1"
    missing_gt = str(c["missingGTinSamplesid"]) == "1"

    ids = AnnotSVIDRegistry()
    g_svlen: Dict[str, object] = {}
    vcf_header = ""
    all_samples: List[str] = []

    print(f"...VCF to BED ({datetime.datetime.now().strftime('%B %d %Y - %H:%M')})")
    to_write: List[str] = []
    written_ids: set = set()           # L_squBrack_SV_ID_Written
    rescue: Dict[str, str] = {}         # TextToWrite_rescue
    rescue_line: Dict[str, int] = {}    # VCFline_rescue
    i = 0
    header_not_present = True
    gt_absent = True
    line_number = 0
    n_unknown_gt = 0
    ls: List[str] = []
    # Variables that keep their value from one VCF line to the next in Tcl.
    chrom2 = pos_vcf2 = end2 = alt_vcf2 = ""
    pos2 = 0
    l2s: List[str] = []

    with _open_vcf(vcf_file) as f:
        for line in _read_lines(f):
            if i == WRITE_CHUNK:
                write_text_in_file("\n".join(to_write), sv_bed)
                to_write = []
                i = 0
            ls = line.split("\t")
            line_number += 1
            if line.startswith("#") or line == "":
                if line.startswith("#CHROM"):
                    for classic in CLASSIC_HEADER:
                        if classic not in ls:
                            raise ConfigExit(f"VCF header file badly formatted:\n{line}\n"
                                             f"({classic} not present)\nExit with error", 2)
                    header_not_present = False
                    if sv_input_info:
                        vcf_header = "SV_type\tSamples_ID\t" + "\t".join(ls[2:])
                    else:
                        vcf_header = "SV_type\tSamples_ID\t" + "\t".join(ls[3:5]) + "\t" + "\t".join(ls[8:])
                    all_samples = ls[9:]
                continue
            i += 1

            if header_not_present:
                raise ConfigExit(f'ERROR:\n{vcf_file}: no VCF header line (prefixed with "#CHROM"). '
                                 "Check your VCF.\nExit with error\n", 2)

            bnd_rescue = False
            tra_rescue = False
            square_bracketed = False
            chrom = _first_sub("chr", "", _idx(ls, 0), re.I)
            pos_vcf = _idx(ls, 1)
            pos = _int(pos_vcf) - 1
            ref = _first_sub("chr", "", _idx(ls, 3))
            alt = _first_sub("chr", "", _idx(ls, 4))
            alt_vcf = alt
            ref_vcf = ref
            ls = _regsub_on_list("CHR2=chr", "CHR2=", ls)
            info = _idx(ls, 7).replace('"', "")

            # Very large DEL written with the sequence in REF: use <DEL>.
            if re.search(r"^[acgtnACGTN]+$", ref, re.I):
                ref_length = len(ref)
                if ref_length > 10000:
                    alt_length = len(alt)
                    if ref[:alt_length] == alt:
                        ref = alt
                        ls = _lset(ls, 3, ref)
                        alt = "<DEL>"
                        ls = _lset(ls, 4, alt)
                        ref_vcf, alt_vcf = ref, alt
                        if not re.search(";SVTYPE=", info) and not re.search("^SVTYPE=", info):
                            info += ";SVTYPE=DEL"
                        if not re.search(";SVLEN=", info) and not re.search("^SVLEN=", info):
                            info += f";SVLEN=-{ref_length - alt_length}"
                        if not re.search(";END=", info) and not re.search("^END=", info):
                            info += f";END={pos + ref_length - alt_length}"
                        ls = _lset(ls, 7, info)

            svlen: object = ""
            end: object = ""
            svtype = cipos = ciend = chr2 = cipos95 = ciend95 = ""
            m = re.search(r"SVLEN=([0-9-]+)", info)
            if m:
                svlen = m.group(1)
            m = re.search(r";END=([0-9]+)", info) or re.search(r"^END=([0-9]+)", info)
            if m:
                end = m.group(1)
            for name in ("SVTYPE", "CIPOS", "CIEND", "CHR2", "CIPOS95", "CIEND95"):
                m = re.search(name + r"=([^;]+)", info)
                if m:
                    value = m.group(1)
                    if name == "SVTYPE":
                        svtype = value
                    elif name == "CIPOS":
                        cipos = value
                    elif name == "CIEND":
                        ciend = value
                    elif name == "CHR2":
                        chr2 = value
                    elif name == "CIPOS95":
                        cipos95 = value
                    else:
                        ciend95 = value

            msq = _SQUARE_BRACKET_ALT.search(alt)
            if alt.startswith("<"):
                # Type2: angle-bracketed notation <...>
                if end == "":
                    end = pos_vcf
                if "<TRA>" in alt:
                    end = pos_vcf
                    svtype = "TRA"
                    svlen = 0
                    if chr2 != "":
                        alt = f"<TRA>_{chr2}_{end}"
                        ls = _lset(ls, 4, alt)
                else:
                    svtype_alt = ""
                    m = re.search(r"^<(.+)>$", alt)
                    if m:
                        svtype_alt = normalize_svtype(m.group(1))
                    svtype = svtype_alt if svtype_alt != "" else normalize_svtype(svtype)
                    if svtype == "INS":
                        m = re.search(r"(INS|INSERTED_SEQUENCE|INSSEQ|SEQ|CONSENSUS)=([acgtnACGTN]+)", info)
                        if m:
                            alt += "_" + m.group(2)
                        elif svlen != "":
                            alt += f"_{svlen}"
                        ls = _lset(ls, 4, alt)

                if svlen == "" and end != "":
                    if svtype == "DEL":
                        svlen = _int(pos_vcf) - _int(end)
                    elif re.search("DUP|INV", svtype):
                        svlen = _int(end) - _int(pos_vcf)
                if svlen != "" and end == "":
                    if svtype == "DEL":
                        end = _int(pos_vcf) - _int(svlen)
                    elif re.search("DUP|INV", svtype):
                        end = _int(svlen) + _int(pos_vcf)

                if svlen != "" and svtype in ("DUP", "DEL", "INS"):
                    if abs(_num(svlen)) < sv_min_size:
                        skip(f"{chrom}_{pos_vcf}_{end}_{svtype}_{ref}_{alt_vcf}: variantLength "
                             f"({abs(_num(svlen))}) < SVminSize ({c['SVminSize']}) (line {line_number})")
                        continue

            elif msq:
                # Type3: square-bracketed notation (one VCF line per breakend)
                base_left, bracket_left, in_chrom, in_start, bracket_right, base_right = msq.groups()
                square_bracketed = True
                if tcl.compare(chrom, in_chrom) != 0:
                    # TRA: both breakends are written in the BED file.
                    ref = "N"
                    end = pos_vcf
                    svtype = "TRA"
                    svlen = 0
                    tra_rescue = True
                    chrom2 = in_chrom
                    pos2 = _int(in_start) - 1
                    pos_vcf2 = in_start
                    end2 = in_start
                    info2 = f"{info};BNDrescue"
                    alt_vcf2 = f"N]{chrom}:{pos_vcf}]"
                    l2 = (f"{chrom2}\t{pos_vcf2}\t{_idx(ls, 2)}\t{ref}\t{alt_vcf2}\t"
                          f"{chr(9).join(ls[5:7])}\t{info2}\t{chr(9).join(ls[8:])}")
                    l2s = [e.replace("\\", "") for e in l2.split("\t")]
                else:
                    if tcl.compare(pos_vcf, in_start) > 0:
                        # Second breakend of a pair: only kept to rescue a missing mate.
                        bnd_rescue = True
                        pos_vcf, in_start = in_start, pos_vcf
                        pos = _int(pos_vcf) - 1
                        cipos, ciend = ciend, cipos
                        info = _first_sub(r"CIPOS=[^;]+", f"CIPOS={cipos}", info)
                        info = _first_sub(r"CIEND=[^;]+", f"CIEND={ciend}", info)
                        cipos95, ciend95 = ciend95, cipos95
                        info = _first_sub(r"CIPOS95=[^;]+", f"CIPOS95={cipos95}", info)
                        info = _first_sub(r"CIEND95=[^;]+", f"CIEND95={ciend95}", info)
                        info += ";BNDrescue"
                        if len(base_left) > 1:
                            base_left = _tcl_regsub_safe("^" + ref, "", base_left, re.I) + "N"
                        elif len(base_right) > 1:
                            base_right = "N" + _tcl_regsub_safe(ref + "$", "", base_right, re.I)
                        elif len(base_left) == 1:
                            base_left = "N"
                        elif len(base_right) == 1:
                            base_right = "N"
                        ref = "N"
                        base_left, base_right = base_right, base_left
                        if bracket_left == "[":
                            bracket_left = bracket_right = "]"
                        else:
                            bracket_left = bracket_right = "["
                        alt = f"{base_left}{bracket_left}{in_chrom}:{in_start}{bracket_right}{base_right}"
                        line = (f"{chrom}\t{pos_vcf}\t{_idx(ls, 2)}\t{ref}\t{alt_vcf}\t"
                                f"{chr(9).join(ls[5:7])}\t{info}\t{chr(9).join(ls[8:])}")
                        ls = line.split("\t")

                    if (len(base_left) > 1 or len(base_right) > 1) and tcl.compare(in_start, str(_int(pos_vcf) + 1)) == 0:
                        # INS
                        if len(base_left) > 1:
                            inserted = _tcl_regsub_safe("^" + ref, "", base_left, re.I)
                        else:
                            inserted = _tcl_regsub_safe(ref + "$", "", base_right, re.I)
                        svlen = len(inserted)
                        svtype = "INS"
                        end = pos_vcf
                        if svlen < sv_min_size:
                            skip(f"{chrom}_{pos_vcf}_{end}_{svtype}_{ref}_{alt_vcf}: variantLength ({svlen}) "
                                 f"< SVminSize ({c['SVminSize']}) (line {line_number})")
                            continue
                        if not re.search("INS|ALL", brack_mode) or svlen > max_size:
                            svlen, svtype, end = 0, "BND", pos_vcf
                        ls = _lset(ls, 4, alt)
                    elif bracket_right == "]" and base_right != "":
                        # DUP
                        svlen = _int(in_start) - _int(pos_vcf)
                        svtype = "DUP"
                        end = in_start
                        if svlen < sv_min_size:
                            skip(f"{chrom}_{pos_vcf}_{end}_{svtype}_{ref}_{alt_vcf}: variantLength ({svlen}) "
                                 f"< SVminSize ({c['SVminSize']}) (line {line_number})")
                            continue
                        if not re.search("DUP|ALL", brack_mode) or svlen > max_size:
                            svlen, svtype, end = 0, "BND", pos_vcf
                    elif (bracket_right == "]" and base_left != "") or (bracket_right == "[" and base_right != ""):
                        # INV
                        svlen = _int(in_start) - _int(pos_vcf)
                        svtype = "INV"
                        end = in_start
                        if not re.search("INV|ALL", brack_mode) or svlen > max_size:
                            svlen, svtype, end = 0, "BND", pos_vcf
                    elif bracket_right == "[" and base_left != "":
                        # DEL (negative svlen)
                        svlen = _int(pos_vcf) - _int(in_start)
                        svtype = "DEL"
                        end = in_start
                        if abs(svlen) < sv_min_size:
                            skip(f"{chrom}_{pos_vcf}_{end}_{svtype}_{ref}_{alt_vcf}: variantLength ({abs(svlen)}) "
                                 f"< SVminSize ({c['SVminSize']}) (line {line_number})")
                            continue
                        if not re.search("DEL|ALL", brack_mode) or svlen > max_size:
                            svlen, svtype, end = 0, "BND", pos_vcf

            elif re.search(r"^[acgtnACGTN.*]+$", ref + alt, re.I):
                # Type1: sequence notation
                refbis = re.sub(r"[*.]", "", ref)
                altbis = re.sub(r"[*.]", "", alt)
                length1: object = len(altbis) - len(refbis)
                if abs(length1) < sv_min_size:
                    skip(f"{chrom}_{pos_vcf}_{end}_{svtype}_{ref}_{alt_vcf}: variantLength ({abs(length1)}) "
                         f"< SVminSize ({c['SVminSize']}) (line {line_number})")
                    continue
                if length1 > 0:
                    if end == "":
                        end = pos_vcf
                    if svtype == "":
                        svtype = "INS"
                else:
                    if end == "":
                        end = pos - length1
                    if svtype == "":
                        svtype = "DEL"
                if "." in ref or "." in alt:
                    # GRIDSS single breakend: length unknown
                    length1 = ""
                elif len(refbis) < 10000 and len(altbis) < 10000:
                    if refbis.lower() not in altbis.lower() and altbis.lower() not in refbis.lower():
                        length1 = ""  # complex SV
                else:
                    length1 = ""
                if svlen == "":
                    svlen = length1
            else:
                skip(f"{chrom}_{pos_vcf}_{end}_{svtype}_{ref}_{alt_vcf}: not a known SV format (line {line_number})")
                continue

            if _lt(end, pos):
                end, pos = pos, end
            if str(end) == str(pos):
                end = _int(pos) + 1

            if include_ci:
                cipos_left = _idx(cipos.split(",") if cipos != "" else [], 0)
                ciend_parts = ciend.split(",") if ciend != "" else []
                ciend_right = ciend_parts[-1] if ciend_parts else ""
                if _lt(cipos_left, 0):
                    pos = _ci_incr(pos, cipos_left)
                    if tcl.compare(str(pos), "0") <= 0:
                        pos = 1
                if tcl.compare(ciend_right, "0") > 0:
                    end = _ci_incr(end, ciend_right)

            if chrom not in ORDERED_CHROMS:
                skip(f'{chrom}_{pos_vcf}_{end}_{svtype}_{ref}_{alt_vcf}: chromosome "{chrom}" unknown (line {line_number})')
                continue
            if end == "":
                skip(f"{chrom}_{pos_vcf}_{end}_{svtype}_{ref}_{alt_vcf}: END of the SV not defined (line {line_number})")
                continue

            if bnd_rescue:
                sv_id = ids.get(f"{chrom}_{pos_vcf}_{end}_{svtype}", ref, alt)
                if sv_id in written_ids:
                    skip(f"{sv_id}: reciprocal breakend (line {line_number})")
                    continue

            if square_bracketed and svtype == "TRA":
                if f"{chrom2}_{pos_vcf2}_{end2}_{svtype}_{ref}_<TRA>" in written_ids:
                    tra_rescue = False

            # Samples_ID
            fmt = _idx(ls, 8)
            fmt_fields = fmt.split(":") if fmt != "" else []
            i_gt = fmt_fields.index("GT") if "GT" in fmt_fields else -1
            if i_gt == -1:
                samples_id = list(all_samples)
            else:
                gt_absent = False
                samples_id = []
                for isample, value in enumerate(ls[9:]):
                    gt = _idx(value.split(":") if value != "" else [], i_gt)
                    if gt not in ("0/0", "0|0"):
                        if gt in ("./.", ".|."):
                            n_unknown_gt += 1
                            if missing_gt:
                                samples_id.append(_idx(all_samples, isample))
                        else:
                            samples_id.append(_idx(all_samples, isample))
            sid = ",".join(samples_id)

            deb = f"{chrom}_{pos_vcf}_{end}_{svtype}"
            annotsv_id = ids.get(deb, ref, alt)
            if sv_input_info:
                if bnd_rescue:
                    rescue[annotsv_id] = (f"{chrom}\t{pos}\t{end}\t{svtype}\t{sid}\t{_idx(ls, 2)}\t"
                                          f"{ref_vcf}\t{alt_vcf}\t" + "\t".join(ls[5:]))
                    rescue_line[annotsv_id] = line_number
                else:
                    to_write.append(f"{chrom}\t{pos}\t{end}\t{svtype}\t{sid}\t" + "\t".join(ls[2:]))
                    if square_bracketed:
                        written_ids.add(annotsv_id)
                if tra_rescue:
                    annotsv_id2 = ids.get(f"{chrom2}_{pos_vcf2}_{end2}_{svtype}", "N", alt_vcf2)
                    rescue[annotsv_id2] = (f"{chrom2}\t{pos2}\t{end2}\t{svtype}\t{sid}\t{_idx(l2s, 2)}\t"
                                           f"{ref_vcf}\t{alt_vcf}\t" + "\t".join(l2s[5:]))
                    rescue_line[annotsv_id2] = line_number
            else:
                if bnd_rescue:
                    rescue[annotsv_id] = (f"{chrom}\t{pos}\t{end}\t{svtype}\t{sid}\t{ref_vcf}\t{alt_vcf}\t"
                                          + "\t".join(ls[8:]))
                    rescue_line[annotsv_id] = line_number
                else:
                    to_write.append(f"{chrom}\t{pos}\t{end}\t{svtype}\t{sid}\t{ref_vcf}\t{alt_vcf}\t"
                                    + "\t".join(ls[8:]))
                    if square_bracketed:
                        if svtype == "TRA":
                            written_ids.add(f"{chrom}_{pos_vcf}_{end}_{svtype}_{ref}_<TRA>")
                        else:
                            written_ids.add(annotsv_id)
                if tra_rescue:
                    annotsv_id2 = ids.get(f"{chrom2}_{pos_vcf2}_{end2}_{svtype}", "N", alt_vcf2)
                    rescue[annotsv_id2] = (f"{chrom2}\t{pos2}\t{end2}\t{svtype}\t{sid}\t{ref_vcf}\t{alt_vcf}\t"
                                           + "\t".join(l2s[8:]))
                    rescue_line[annotsv_id2] = line_number

            if svlen != "":
                g_svlen[annotsv_id] = svlen
                if tra_rescue:
                    g_svlen[annotsv_id2] = 0

    if missing_gt:
        if n_unknown_gt > 1:
            print(f"\t...WARNING: {n_unknown_gt} missing alleles (./. or .|.) found in the GT field of reported SV")
            print('\t            Corresponding samples have been reported in the "Samples_ID" output field\n')
        elif n_unknown_gt == 1:
            print("\t...WARNING: 1 missing allele (./. or .|.) found in the GT field of reported SV")
            print('\t            Corresponding samples have been reported in the "Samples_ID" output field\n')
        else:
            print("\n")

    write_text_in_file("\n".join(to_write), sv_bed)
    to_write = []
    for sv_id, text in rescue.items():
        if sv_id not in written_ids:
            to_write.append(text)
        else:
            skip(f"{sv_id}: reciprocal breakend (line {rescue_line[sv_id]})")
    write_text_in_file("\n".join(to_write), sv_bed)

    if c.get("candidateSnvIndelFiles") and gt_absent:
        print(f"...parsing of {vcf_file}")
        print("\t...WARNING: The SV genotype is not indicated in the FORMAT column under the “GT” field")
        print("\t            => Compound heterozygosity analysis won't be processed!\n")
        c["candidateSnvIndelFiles"] = []
        c["candidateSnvIndelSamples"] = []

    if vcf_header == "":
        # No "#CHROM" line (only possible without any data line).
        if sv_input_info:
            vcf_header = "SV_type\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join([""] * len(ls[9:]))
        else:
            vcf_header = "SV_type\tREF\tALT\t" + "\t".join([""] * len(ls[8:]))

    if is_an_empty_file(sv_bed):
        raise ConfigExit("############################################################################\n"
                         "No SV to annotate in the SVinputFile - Exit without error.\n"
                         "############################################################################", 0)

    return VCFtoBEDResult(bed_file=sv_bed, vcf_header=vcf_header, svlen=g_svlen, unannotated_file=unannotated)


def _lset(items: List[str], index: int, value: str) -> List[str]:
    """``lset`` (Tcl errors when the index is past the end; we pad instead)."""
    items = list(items)
    while len(items) <= index:
        items.append("")
    items[index] = value
    return items
