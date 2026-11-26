"""Output writing: first part of the port of ``AnnotSV-write.tcl``.

``organize_annotation`` is the walking skeleton of ``OrganizeAnnotation``:
it reads the "full and split" file made by ``genes.genes_annotation`` and
writes the output TSV with the core columns, from ``AnnotSV_ID`` to
``RE_gene``:

    AnnotSV_ID SV_chrom SV_start SV_end SV_length [input columns]
    Annotation_mode Gene_name Closest_left Closest_right Gene_count Tx
    Tx_version Tx_start Tx_end Overlapped_tx_length Overlapped_CDS_length
    Overlapped_CDS_percent Frameshift Exon_count Location Location2
    Dist_nearest_SS Nearest_SS_type Intersect_start Intersect_end RE_gene

The gene-based columns (``gene_based.GeneBasedAnnotations``) follow
``RE_gene``, and ``-candidateGenesFiltering`` is applied. The other column
groups (cytoband, pathogenic/benign SV, breakpoints, ranking...) are ported
in later modules and go between ``RE_gene`` and the gene-based columns. Closest genes
and regulatory elements come in through the ``closest_genes`` and
``regulatory`` arguments and stay empty until those modules exist.
"""

from __future__ import annotations

import glob
import itertools
import os
import re
from typing import Callable, Dict, List, Optional

from . import tcl
from .general import (
    AnnotSVIDRegistry,
    check_bed,
    first_line_from_file,
    normalize_svtype,
    replace_text_in_file,
    write_text_in_file,
)
from .filtered_vcf import FilteredVCF
from .gene_based import is_candidate
from .genes import GenesResult, TranscriptVersions, unix_sort_bed
from .snv_indel import SnvIndelCounts

GENE_COLUMNS = ("Gene_name\tClosest_left\tClosest_right\tGene_count\tTx\tTx_version\tTx_start\tTx_end"
                "\tOverlapped_tx_length\tOverlapped_CDS_length\tOverlapped_CDS_percent\tFrameshift"
                "\tExon_count\tLocation\tLocation2\tDist_nearest_SS\tNearest_SS_type\tIntersect_start"
                "\tIntersect_end\tRE_gene")

_COORD_RE = re.compile(r"[0-9XYMT]+:[0-9]+-[0-9]+")


# --------------------------------------------------------------------------
# BED -> VCF coordinates
# --------------------------------------------------------------------------
def switch_coordinates_from_bed_to_vcf(bed_coord: str) -> str:
    """switchCoordinatesFromBEDtoVCF: ``chr1:12356-12359`` -> ``chr1:12357-12359``."""
    parts = tcl.split_chars(bed_coord, ":|-")
    if len(parts) != 3:
        return bed_coord
    return f"{parts[0]}:{int(parts[1]) + 1}-{parts[2]}"


def switch_all_coordinates_from_bed_to_vcf_in_line(line: str) -> str:
    """switchAllCoordinatesFromBEDtoVCFinLine: +1 on the start in AnnotSV_ID,
    SV_start and every ``chrom:start-end`` of the line."""
    sv_start_bed = re.match(r"[^\t]+\t[^\t]+\t([0-9]+)", line).group(1)
    sv_start_vcf = str(int(sv_start_bed) + 1)

    annotsv_id = tcl.lindex(line, 0)
    parts = annotsv_id.split("_")
    get = lambda k: parts[k] if k < len(parts) else ""  # noqa: E731
    updated = f"{get(0)}_{sv_start_vcf}_{get(2)}_{get(3)}_{get(4)}"
    line = updated + line[len(annotsv_id):]

    chrom = tcl.lindex(line, 1)
    m = re.compile(sv_start_bed).search(line, len(updated) + len(chrom) + 1)
    if m:
        line = line[:m.start()] + sv_start_vcf + line[m.end():]

    i = 0
    while True:
        m = _COORD_RE.search(line, i)
        if not m:
            break
        line = line[:m.start()] + switch_coordinates_from_bed_to_vcf(m.group(0)) + line[m.end():]
        i = m.end()
    return line


# --------------------------------------------------------------------------
# Gene columns of one line
# --------------------------------------------------------------------------
def _lt(a, b) -> bool:
    return tcl.compare(str(a), str(b)) < 0


def _le(a, b) -> bool:
    return tcl.compare(str(a), str(b)) <= 0


def _is_integer(value) -> bool:
    """Tcl ``string is integer`` (the empty string is accepted)."""
    return isinstance(value, int) or value == "" or re.fullmatch(r"\s*[+-]?[0-9]+\s*", str(value)) is not None


def _frameshift(cds_l: str, svtype_norm: str, svlen) -> str:
    if not _is_integer(cds_l) or cds_l == "":
        return "NA"
    cds_l = int(cds_l)
    if svtype_norm in ("DEL", "DUP"):
        if cds_l == 0:
            return "NA"           # outside the CDS
        return "no" if cds_l % 3 == 0 else "yes"
    if svtype_norm == "INS":
        if svlen is None or not _is_integer(svlen) or svlen == "":
            return "NA"           # insertion length unknown
        return "NA" if int(svlen) % 3 == 0 else "yes"
    return "NA"                   # TRA, INV and other SV types


def _ss_in_exon(d_a: int, d_b: int, i: int, nb_ex: int, plus: bool):
    """Nearest splice site when the breakpoint is in exon ``i``."""
    if i == 1:                    # only 1 splice site on the first exon
        return d_b, "5'" if plus else "3'"
    if i == nb_ex:                # only 1 splice site on the last exon
        return d_a, "3'" if plus else "5'"
    if d_a < d_b:
        return d_a, "3'" if plus else "5'"
    return d_b, "5'" if plus else "3'"


def _ss_in_intron(d_b: int, d_a: int, plus: bool):
    if d_b < d_a:
        return d_b, "5'" if plus else "3'"
    return d_a, "3'" if plus else "5'"


def split_line_gene_columns(ls: List[str], sv_left: str, sv_right: str, svtype: str, svlen,
                            tx_version: Callable[[str], str]) -> Dict[str, str]:
    """Gene columns of a "split" line (SV line + gene line + lengths + "split")."""
    strand = ls[-11]
    exon_starts, exon_ends = ls[-6], ls[-5]
    cds_start, cds_end = ls[-8], ls[-7]
    cds_l = ls[-4]
    plus = strand == "+"
    col = {
        "Gene_name": ls[-10], "Gene_count": "", "Tx": ls[-9], "Tx_version": tx_version(ls[-9]),
        "Tx_start": ls[-13], "Tx_end": ls[-12], "Overlapped_tx_length": ls[-2],
        "Overlapped_CDS_length": cds_l, "Overlapped_CDS_percent": ls[-3],
    }
    starts, ends = exon_starts.split(","), exon_ends.split(",")
    col["Exon_count"] = str(len(starts) - 1)
    svtype_norm = normalize_svtype(svtype)
    frameshift = _frameshift(cds_l, svtype_norm, svlen)

    # Location and nearest splice site
    nb_ex = len(starts) - 1
    tx_left = starts[0]
    tx_right = ends[-2] if len(ends) >= 2 else ""
    loc_start = loc_end = ""
    dist_left = dist_right = ""
    type_left = type_right = ""
    if _le(sv_left, tx_left):     # SV begins before tx start (or tx end for strand "-")
        if plus:
            loc_start = "txStart"
        else:
            loc_end = "txEnd"
    i = 0
    previous_b = ends[0]
    left, right = int(sv_left), int(sv_right)
    for a, b in itertools.zip_longest(starts, ends, fillvalue=""):
        if a == "" or b == "":
            continue
        i += 1
        a, b = int(a), int(b)
        # SV left
        if left < b:
            if left > a:          # in an exon
                dist_left, type_left = _ss_in_exon(left - a, b - left, i, nb_ex, plus)
                if plus:
                    loc_start = f"exon{i}"
                else:
                    loc_end = f"exon{nb_ex - i + 1}"
            elif left >= int(previous_b):   # in an intron
                dist_left, type_left = _ss_in_intron(left - int(previous_b), a - left, plus)
                if plus:
                    loc_start = f"intron{i - 1}"
                else:
                    loc_end = f"intron{nb_ex - i + 1}"
        # SV right
        if right < b:
            if right > a:         # in an exon
                dist_right, type_right = _ss_in_exon(right - a, b - right, i, nb_ex, plus)
                if plus:
                    loc_end = f"exon{i}"
                else:
                    loc_start = f"exon{nb_ex - i + 1}"
                break
            elif right >= int(previous_b):  # in an intron
                dist_right, type_right = _ss_in_intron(right - int(previous_b), a - right, plus)
                if plus:
                    loc_end = f"intron{i - 1}"
                else:
                    loc_start = f"intron{nb_ex - i + 1}"
                break
        previous_b = b
    if loc_end == "" or loc_start == "":   # SV finishes after tx end
        if plus:
            loc_end = "txEnd"
        else:
            loc_start = "txStart"
    location = f"{loc_start}-{loc_end}"
    dist, ss_type = "", ""
    if location != "txStart-txEnd":
        if dist_right == "":
            if dist_left != "":
                dist, ss_type = dist_left, type_left
        elif dist_left == "":
            dist, ss_type = dist_right, type_right
        elif dist_right < dist_left:
            dist, ss_type = dist_right, type_right
        else:
            dist, ss_type = dist_left, type_left

    # Location2 (UTR / CDS). Some genes have CDSstart=CDSend=txEnd (no CDS)
    if cds_start == cds_end:
        region_left = region_right = "UTR"
    else:
        def region(pos):
            if _lt(pos, cds_start):
                return "5'UTR" if plus else "3'UTR"
            if _lt(pos, cds_end):
                return "CDS"
            return "3'UTR" if plus else "5'UTR"
        region_left, region_right = region(sv_left), region(sv_right)
    region_start, region_end = (region_left, region_right) if plus else (region_right, region_left)
    location2 = f"{region_start}-{region_end}"
    m = re.search(r"([^-]+)-([^-]+)", location2)
    if m and m.group(1) == m.group(2):
        location2 = m.group(1)
    # Insertion outside the CDS (as in Tcl, Location2 never holds "exon")
    if svtype_norm == "INS" and "exon" not in location2:
        frameshift = "NA"

    col.update({
        "Frameshift": frameshift, "Location": location, "Location2": location2,
        "Dist_nearest_SS": str(dist), "Nearest_SS_type": ss_type,
        "Intersect_start": tx_left if _lt(sv_left, tx_left) else sv_left,
        "Intersect_end": sv_right if _lt(sv_right, tx_right) else tx_right,
    })
    return col


def full_line_gene_columns(sv_key: str, genes_by_sv: Dict[str, str]) -> Dict[str, str]:
    gene_name = genes_by_sv.get(sv_key, "")
    count = str(len(gene_name.split(";"))) if sv_key in genes_by_sv else "0"
    col = {k: "" for k in GENE_COLUMNS.split("\t")}
    col["Gene_name"], col["Gene_count"] = gene_name, count
    return col


# --------------------------------------------------------------------------
# OrganizeAnnotation (core columns)
# --------------------------------------------------------------------------
def _header_svinput_file(bed_file: str) -> str:
    return re.sub(r"(.formatted)?(.sorted)?.bed$", lambda _m: ".header.tsv", bed_file, count=1, flags=re.I)


def build_header(config, vcf_header: Optional[str], cytoband: str = "") -> (str, int):
    """The header up to ``RE_gene`` and ``theBEDlength``. As in Tcl, updates
    ``SVinputInfo`` / ``formatColNumber`` (VCF) or ``samplesidTSVcol`` (BED)."""
    c = config
    header = "AnnotSV_ID\tSV_chrom\tSV_start\tSV_end\tSV_length"
    bed_length = 0
    if vcf_header is not None:
        header += f"\t{vcf_header}"
        bed_length = len(header.split("\t")) - 2
        c["formatColNumber"] = 10 if str(c["SVinputInfo"]) == "1" else 6
        if "Samples_ID" in c.get("outputColHeader", []):
            c["formatColNumber"] += 1
        c["SVinputInfo"] = "1"
    elif str(c["SVinputInfo"]) == "1":
        header_file = _header_svinput_file(c["bedFile"])
        bed_length = len(first_line_from_file(c["bedFile"]).split("\t"))
        from_user = first_line_from_file(header_file).split("\t")
        if len(from_user) != bed_length:
            print(f"Numbers of columns from {c['SVinputFile']} and {header_file} are different "
                  f"({bed_length} != {len(from_user)})")
            print(f"=> Can not report: {tcl.format_list(from_user)}")
        else:
            if int(c["svtBEDcol"]) != -1:
                from_user[int(c["svtBEDcol"])] = "SV_type"
            if int(c["samplesidBEDcol"]) != -1:
                from_user[int(c["samplesidBEDcol"])] = "Samples_ID"
            k = 1
            for name in from_user[3:]:
                if name != "":
                    header += f"\t{name}"
                else:
                    header += f"\tuser#{k}"
                    k += 1
    else:
        if int(c.get("svtTSVcol", -1)) != -1:
            header += "\tSV_type"
        if int(c.get("samplesidTSVcol", -1)) == -1:
            names = first_line_from_file(_header_svinput_file(c["bedFile"])).split("\t")
            c["samplesidTSVcol"] = str(names.index("Samples_ID") if "Samples_ID" in names else -1)
        header += "\tSamples_ID"
    header += "\tAnnotation_mode" + cytoband + "\t" + GENE_COLUMNS
    return header, bed_length


def write_users_bed(config, full_and_split_file: str) -> str:
    """The sorted "chrom start end" BED of every full SV and every SV/gene
    intersection (``g_AnnotSV(fullAndSplitBedFile)``), input of the region
    annotations."""
    c = config
    users = f"{c['outputDir']}/{os.path.basename(c['bedFile'])}.users.bed"
    if os.path.exists(users):
        os.remove(users)
    text = set()
    with open(full_and_split_file, "r", newline="") as f:
        for raw in f:
            line = raw[:-1] if raw.endswith("\n") else raw
            if line == "":
                continue
            ls = line.split("\t")
            if ls[-1] == "split":
                tx_left = ls[-6].split(",")[0]
                ends = ls[-5].split(",")
                tx_right = ends[-2] if len(ends) >= 2 else ""
                start = tx_left if _lt(ls[1], tx_left) else ls[1]
                end = ls[2] if _lt(ls[2], tx_right) else tx_right
                text.add(f"{ls[0]}\t{start}\t{end}")
            else:
                text.add("\t".join(ls[0:3]))
    write_text_in_file("\n".join(sorted(text)), users)
    check_bed(users, c["outputDir"])
    formatted = re.sub(r".bed$", lambda _m: ".formatted.bed", users, count=1, flags=re.I)
    sorted_users = f"{c['outputDir']}/{os.path.basename(c['bedFile'])}.users.sorted.bed"
    if os.path.exists(sorted_users):
        os.remove(sorted_users)
    unix_sort_bed(formatted, sorted_users, "OrganizeAnnotation")
    os.remove(formatted)
    c["fullAndSplitBedFile"] = sorted_users
    return sorted_users


def organize_annotation(config, genes: GenesResult, vcf_header: Optional[str] = None,
                        svlen: Optional[Dict[str, object]] = None,
                        closest_genes: Optional[Callable[[str, str, str], List[str]]] = None,
                        regulatory: Optional[Dict[str, str]] = None,
                        header_file_to_remove: bool = False,
                        gene_based=None, regions=None, pathogenic=None, benign=None) -> str:
    """OrganizeAnnotation, core columns only. Returns the output file path.

    ``vcf_header`` and ``svlen`` come from the SV input step (``VCFheader``,
    ``g_SVLEN``). As in Tcl, the AnnotSV_IDs are numbered afresh here (VCFsToBED
    unsets ``g_ID`` / ``g_deb`` when it is done).
    ``closest_genes(chrom, pos, "left"|"right")`` stands for
    closestGenesAnnotation and ``regulatory`` for ``g_re``. ``gene_based`` is a
    ``gene_based.GeneBasedAnnotations`` (the "Gene-based" columns, which come
    right after ``RE_gene`` until the column groups between them are ported).
    ``regions`` is a ``regions.RegionAnnotations``: CytoBand after
    ``Annotation_mode``, then after ``RE_gene`` the user SVincludedInFt BED
    files, TAD, the user FtIncludedInSV, AnyOverlap and BNDproximity BED files
    and the breakpoint columns (GC content, Repeat, Gap, SegDup, ENCODE
    blacklist), before the gene-based columns. ``pathogenic`` is a
    ``pathogenic.PathogenicAnnotations``: the pathogenic SV, partially
    overlapped pathogenic SV and pathogenic SNV/indel columns, right after
    ``RE_gene``. ``benign`` is a ``benign.BenignAnnotations``: the benign SV
    and partially overlapped benign SV columns, right after the pathogenic
    ones. With ``snvIndelFiles`` / ``candidateSnvIndelFiles`` set, the
    SNV/indel counts (``snv_indel.SnvIndelCounts``) and the compound_htz
    columns (``config.filtered_vcf``) come right after TAD.
    """
    c = config
    svlen = svlen if svlen is not None else {}
    ids = AnnotSVIDRegistry()
    closest_genes = closest_genes or (lambda chrom, pos, side: [])
    regulatory = regulatory or {}
    tx_version = TranscriptVersions(c)
    output_file = f"{c['outputDir']}/{c['outputFile']}"

    REGION_KEYS = ("SVincludedInFt", "TAD", "FtIncludedInSV", "AnyOverlap", "BNDproximity", "breakpoints")
    cyto_header = regions.cytoband_header() if regions is not None else ""
    header, bed_length = build_header(c, vcf_header, cyto_header)
    if pathogenic is not None:
        header += pathogenic.header()
    if benign is not None:
        header += benign.header()
    # SNV/indel counts (VCFannotation) and compound htz (filteredVCFannotation):
    # after TAD, before the user FtIncludedInSV BED files.
    snv_counts = SnvIndelCounts(c) if c.get("snvIndelFiles") else None
    compound = None
    snv_header = snv_counts.header() if snv_counts is not None else ""
    if c.get("candidateSnvIndelFiles"):
        compound = c.filtered_vcf if c.filtered_vcf is not None else FilteredVCF(c)
        c.filtered_vcf = compound
        snv_header += "".join(f"\tcompound_htz({s})" for s in c["candidateSnvIndelSamples"])
    if regions is not None:
        parts = regions.header_parts()
        header += "".join(parts[k] + (snv_header if k == "TAD" else "") for k in REGION_KEYS)
    else:
        header += snv_header
    if gene_based is not None:
        header += gene_based.header()
    replace_text_in_file(header, output_file)
    hcols = header.split("\t")
    i_ref = (hcols.index("REF") if "REF" in hcols else -1) - 2
    i_alt = (hcols.index("ALT") if "ALT" in hcols else -1) - 2

    write_users_bed(c, genes.full_and_split_file)

    svt_col, sid_col = int(c["svtBEDcol"]), int(c["samplesidBEDcol"])
    svinput_info = str(c["SVinputInfo"]) == "1"
    order: List[str] = []
    text_by_id: Dict[str, List[str]] = {}
    annotsv_id = ""

    with open(genes.full_and_split_file, "r", newline="") as f:
        for raw in f:
            line = raw[:-1] if raw.endswith("\n") else raw
            if line == "":
                continue
            ls = line.split("\t")
            chrom, sv_left, sv_right = ls[0], ls[1], ls[2]
            mode = ls[-1]
            svtype = ls[svt_col] if svt_col != -1 and svt_col < len(ls) else ""
            samplesid = ls[sid_col] if sid_col != -1 and sid_col < len(ls) else ""

            if mode == "split":
                # Tcl reads g_SVLEN with the AnnotSV_ID of the previous line (the "full" line of this SV)
                col = split_line_gene_columns(ls, sv_left, sv_right, svtype, svlen.get(annotsv_id),
                                              tx_version)
                col["Closest_left"] = col["Closest_right"] = ""
                col["RE_gene"] = ""
            else:
                col = full_line_gene_columns("\t".join(ls[0:3]), genes.genes_by_sv)
                overlapped = col["Gene_name"].split(";")
                left = closest_genes(chrom, sv_left, "left")
                col["Closest_left"] = next((g for g in reversed(left) if g not in overlapped), "")
                right = closest_genes(chrom, sv_right, "right")
                col["Closest_right"] = next((g for g in right if g not in overlapped), "")
                col["RE_gene"] = regulatory.get(f"{chrom}\t{sv_left}\t{sv_right}", "")

            # AnnotSV_ID and SV length
            sv_start_bed = sv_left
            sv_start_vcf = int(sv_start_bed) + 1
            ref = (ls[i_ref] if 0 <= i_ref < len(ls) else "").replace("chr", "", 1)
            alt = (ls[i_alt] if 0 <= i_alt < len(ls) else "").replace("chr", "", 1)
            svtype_tmp = re.sub(r"[ .():;]", "_", svtype)
            annotsv_id = ids.get(f"{chrom}_{sv_start_vcf}_{sv_right}_{svtype_tmp}", ref, alt)
            if annotsv_id in svlen:
                sv_length = svlen[annotsv_id]
            else:
                norm = normalize_svtype(svtype)
                if "DEL" in norm:
                    sv_length = int(sv_start_bed) - int(sv_right)
                elif re.search("DUP|INV", norm):
                    sv_length = int(sv_right) - int(sv_start_bed)
                elif "TRA" in norm:
                    sv_length = 0
                    svlen[annotsv_id] = 0
                else:
                    sv_length = ""

            if svinput_info:
                to_add = ls[0:bed_length]
                to_add.insert(3, str(sv_length))
                text = f"{annotsv_id}\t" + "\t".join(to_add)
            else:
                text = f"{annotsv_id}\t{chrom}\t{sv_left}\t{sv_right}\t{sv_length}"
                if svt_col != -1:
                    text += f"\t{svtype}"
                if sid_col != -1:
                    text += f"\t{samplesid}"
            if mode == "split":
                r_start, r_end = col["Intersect_start"], col["Intersect_end"]
            else:
                r_start, r_end = sv_left, sv_right
            p_text = pathogenic.text(chrom, r_start, r_end, sv_left, sv_right, mode) if pathogenic is not None else ""
            if benign is not None:
                p_text += benign.text(chrom, r_start, r_end, sv_left, sv_right, mode, col["Gene_name"])
            snv_text = ""
            if snv_counts is not None:
                snv_text += "\t" + snv_counts.text(chrom, r_start, r_end, svtype)
            if compound is not None:
                if mode == "split":
                    snv_text += "\t" + compound.annotate(chrom, col["Tx_start"], col["Tx_end"], ls, header)
                else:
                    snv_text += "\t" + compound.annotate("FULL")
            if regions is not None:
                cyto = regions.cytoband_text(chrom, r_start, r_end)
                text += f"\t{mode}" + (f"\t{cyto}" if cyto is not None else "")
                text += "\t" + "\t".join(col[k] for k in GENE_COLUMNS.split("\t")) + p_text
                r_texts = regions.texts(chrom, r_start, r_end, sv_left, sv_right, mode)
                text += "".join(r_texts[k] + (snv_text if k == "TAD" else "") for k in REGION_KEYS)
            else:
                text += f"\t{mode}\t" + "\t".join(col[k] for k in GENE_COLUMNS.split("\t")) + p_text + snv_text
            if gene_based is not None and gene_based.enabled:
                text += "\t" + "\t".join(gene_based.text(col["Gene_name"], mode))

            if annotsv_id not in text_by_id:
                order.append(annotsv_id)
                text_by_id[annotsv_id] = []
            if not (c["annotationMode"] == "full" and mode == "split"):
                text_by_id[annotsv_id].append(text)

    i_mode = hcols.index("Annotation_mode")
    i_gene = hcols.index("Gene_name")
    candidate_filtering = str(c.values.get("candidateGenesFiltering", "0")).lower() in ("1", "yes", "true", "on")
    lines: List[str] = []
    for annotsv_id in order:
        for text in text_by_id[annotsv_id]:
            mode = _col(text, i_mode)
            if c["annotationMode"] == "full" and mode == "split":
                continue
            if c["annotationMode"] == "split" and mode == "full":
                continue
            # To select only the SV annotations overlapping a gene from the "candidateGenesFile"
            if candidate_filtering:
                gene_name = _col(text, i_gene)
                if gene_name == "" or not any(is_candidate(c, g) for g in gene_name.split(";")):
                    continue
            lines.append(switch_all_coordinates_from_bed_to_vcf_in_line(text))
    write_text_in_file("\n".join(lines), output_file)

    # Temporary files
    bed = c["bedFile"]
    to_delete = [genes.full_and_split_file, bed, c["fullAndSplitBedFile"],
                 re.sub(r"(.formatted)?.bed$", lambda _m: ".breakpoints.bed", bed, count=1, flags=re.I)]
    to_delete += glob.glob(f"{c['outputDir']}/*_AnnotSV_inputSVfile.bed")
    if header_file_to_remove:
        to_delete.append(_header_svinput_file(bed))
    if c.get("NAbedFile"):
        to_delete.append(c["NAbedFile"])
    for svtype in ("Gain", "Loss"):
        tmp = re.sub(r"(.formatted)?.bed$", lambda _m: f".intersect.po_benign-{svtype}", bed, count=1, flags=re.I)
        to_delete.append(f"{c['outputDir']}/{os.path.basename(tmp)}")
    for path in to_delete:
        if os.path.isfile(path):
            os.remove(path)

    print(f"\n...output columns annotation:")
    print("\t" + header.replace("\t", ";") + "\n")
    return output_file


def _col(text: str, i: int) -> str:
    ls = text.split("\t")
    return ls[i] if i < len(ls) else ""
