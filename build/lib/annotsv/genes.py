"""Gene annotation: port of ``AnnotSV-genes.tcl``.

* ``check_genes_refseq_file`` / ``check_genes_ensembl_file`` prepare the
  ``genes.<tx>.sorted.bed`` file (``checkGenesRefSeqFile`` /
  ``checkGenesENSEMBLfile``).
* ``TranscriptVersions`` answers ``transcriptVersionAnnotation``.
* ``genes_annotation`` intersects the SV BED with the genes file and writes
  the "full and split" file read by the output writer (``genesAnnotation``).

The Tcl globals ``g_Lgenes`` (genes overlapped by each SV) and the list of
all overlapped genes (passed to ``regulatoryElementsAnnotation``) are
returned in a ``GenesResult`` instead.
"""

from __future__ import annotations

import datetime
import glob
import itertools
import os
import re
import subprocess
from dataclasses import dataclass, field
from typing import Dict, List, Optional

from . import tcl
from .config import ConfigExit
from .general import (
    check_samplesid_bed_col,
    check_svt_bed_col,
    lines_from_file,
    lines_from_gz_file,
    remove_redundancy_without_sorting,
    sort_on_element,
    write_text_in_file,
)

KARYOTYPE = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15",
             "16", "17", "18", "19", "20", "21", "22", "X", "Y", "M", "MT"]

_TX_VERSION_RE = re.compile(r"([^.]+)\.([0-9_]+)$")


def genes_dir(config) -> str:
    c = config
    return f"{c['annotationsDir']}/Annotations_{c['organism']}/Genes/{c['genomeBuild']}"


def unix_sort_bed(src: str, dst: str, caller: str) -> None:
    """``LC_ALL=C sort -k1,1 -k2,2n src >> dst`` (the Tcl code shells out too)."""
    env = dict(os.environ, LC_ALL="C")
    with open(dst, "ab") as out:
        r = subprocess.run(["sort", "-k1,1", "-k2,2n", src], stdout=out, stderr=subprocess.PIPE, env=env)
    if r.returncode:
        raise ConfigExit(f"-- {caller} --\nsort -k1,1 -k2,2n {src} > {dst}\n"
                         f"{r.stderr.decode(errors='replace')}\nExit with error", 2)


def _split_tx_version(name: str):
    m = _TX_VERSION_RE.search(name)
    if m:
        return m.group(1), m.group(2)
    return name, ""


def _get(ls: List[str], i: int) -> str:
    return ls[i] if -len(ls) <= i < len(ls) else ""


def _missing(path_text: str) -> ConfigExit:
    bar = "#" * 76
    return ConfigExit(f"{bar}\n{path_text}\nPlease check your install - Exit with error.\n{bar}", 2)


# --------------------------------------------------------------------------
# Genes files
# --------------------------------------------------------------------------
def check_genes_refseq_file(config, now: Optional[datetime.datetime] = None) -> str:
    """checkGenesRefSeqFile: create ``genes.RefSeq.sorted.bed`` (and
    ``transcript_version.RefSeq.tsv``) from the UCSC download if needed.
    Sets and returns ``config["genesFile"]``."""
    c = config
    gdir = genes_dir(c)
    build = c["genomeBuild"]
    downloaded: List[str] = []
    if re.search("GRCh37|GRCh38|mm39", build):
        downloaded = glob.glob(f"{gdir}/ncbiRefSeq.txt.gz")
    elif re.search("mm9|mm10", build):
        downloaded = glob.glob(f"{gdir}/refGene.txt.gz")
    elif re.search("CHM13", build):
        downloaded = sorted(glob.glob(f"{gdir}/hs1_curGene*.txt.gz"))
    formatted = f"{gdir}/genes.RefSeq.sorted.bed"
    tx_version_file = f"{gdir}/transcript_version.RefSeq.tsv"

    if not downloaded and not os.path.exists(formatted):
        raise _missing(f'"{c["annotationsDir"]}/Annotations_{c["organism"]}/Genes/{build}/*.txt.gz" '
                       "file doesn't exist")

    if not os.path.exists(formatted):
        prom_dir = (f"{c['annotationsDir']}/Annotations_{c['organism']}/FtIncludedInSV/"
                    f"RegulatoryElements/{build}")
        for prom in glob.glob(f"{prom_dir}/promoter_*bp_RefSeq_{build}.sorted.bed"):
            os.remove(prom)

        stamp = (now or datetime.datetime.now()).strftime("%B %d %Y - %H:%M")
        print(f"\t...creation of {formatted} ({stamp})")
        print("\t   (done only once during the first annotation)")

        by_chrom: Dict[str, List[str]] = {}
        for line in lines_from_gz_file(downloaded[0]):
            by_chrom.setdefault(_get(line.split("\t"), 2), []).append(line)
        sorted_lines: List[str] = []
        for val in KARYOTYPE:
            if f"chr{val}" in by_chrom:
                sorted_lines += sort_on_element(sort_on_element(by_chrom[f"chr{val}"], 5), 4)

        if os.path.exists(tx_version_file):
            os.remove(tx_version_file)
        tx_versions: List[str] = []
        formatted_lines: List[str] = []
        seen = set()
        for line in sorted_lines:
            ls = line.split("\t")
            chrom = _get(ls, 2).replace("chr", "", 1)
            tx_id, tx_version = _split_tx_version(_get(ls, 1))
            tx_versions.append(f"{tx_id}\t{tx_version}")
            out = "\t".join([chrom, _get(ls, 4), _get(ls, 5), _get(ls, 3), _get(ls, 12), tx_id,
                             _get(ls, 6), _get(ls, 7), _get(ls, 9), _get(ls, 10)])
            if out not in seen:
                formatted_lines.append(out)
                seen.add(out)
        for f in downloaded:
            os.remove(f)
        write_text_in_file("\n".join(formatted_lines), f"{formatted}.tmp")
        write_text_in_file("\n".join(tx_versions), tx_version_file)
        unix_sort_bed(f"{formatted}.tmp", formatted, "checkGenesRefSeqFile")
        os.remove(f"{formatted}.tmp")

    c["genesFile"] = formatted
    return formatted


def check_genes_ensembl_file(config, now: Optional[datetime.datetime] = None) -> str:
    """checkGenesENSEMBLfile: strip the transcript versions from
    ``genes.ENSEMBL.sorted.tmp.bed`` into ``genes.ENSEMBL.sorted.bed`` and
    ``transcript_version.ENSEMBL.tsv``. Sets and returns ``config["genesFile"]``."""
    c = config
    gdir = genes_dir(c)
    tmp = f"{gdir}/genes.ENSEMBL.sorted.tmp.bed"
    formatted = f"{gdir}/genes.ENSEMBL.sorted.bed"
    tx_version_file = f"{gdir}/transcript_version.ENSEMBL.tsv"

    if not os.path.exists(formatted) and not os.path.exists(tmp):
        raise _missing(f'"{formatted}" doesn\'t exist')

    if os.path.exists(tmp) and not os.path.exists(formatted):
        stamp = (now or datetime.datetime.now()).strftime("%B %d %Y - %H:%M")
        print(f"\t...creation of {formatted} ({stamp})")
        print("\t   (done only once during the first annotation)")
        tx_versions, to_write = [], []
        for line in lines_from_file(tmp):
            ls = line.split("\t")
            tx_id, tx_version = _split_tx_version(_get(ls, 5))
            tx_versions.append(f"{tx_id}\t{tx_version}")
            to_write.append("\t".join(ls[0:5] + [tx_id] + ls[6:]))
        write_text_in_file("\n".join(tx_versions), tx_version_file)
        write_text_in_file("\n".join(to_write), formatted)
        os.remove(tmp)

    c["genesFile"] = formatted
    return formatted


def check_genes_file(config) -> str:
    """The ``-tx`` switch of ``bin/AnnotSV``."""
    if config["tx"] == "RefSeq":
        return check_genes_refseq_file(config)
    return check_genes_ensembl_file(config)


class TranscriptVersions:
    """transcriptVersionAnnotation (``g_transcriptVersion`` loaded lazily)."""

    def __init__(self, config) -> None:
        self.config = config
        self._versions: Optional[Dict[str, str]] = None if config["organism"] == "Human" else {}

    def _load(self) -> Dict[str, str]:
        c = self.config
        path = f"{genes_dir(c)}/transcript_version.{c['tx']}.tsv"
        if not os.path.exists(path):
            raise _missing(f'"{path}" doesn\'t exist')
        versions: Dict[str, str] = {}
        for line in lines_from_file(path):
            ls = line.split("\t")
            versions[_get(ls, 0)] = _get(ls, 1)
        return versions

    def __call__(self, transcript_id: str) -> str:
        if self._versions is None:
            self._versions = self._load()
        return self._versions.get(transcript_id, "")


# --------------------------------------------------------------------------
# genesAnnotation
# --------------------------------------------------------------------------
@dataclass
class GenesResult:
    full_and_split_file: str
    # g_Lgenes: "chrom\tstart\tend" -> "GENE1;GENE2" (genes overlapped, coordinate order)
    genes_by_sv: Dict[str, str] = field(default_factory=dict)
    # lsort -unique of every gene overlapped (input of regulatoryElementsAnnotation)
    all_genes_overlapped: List[str] = field(default_factory=list)


def _cds_lengths(exon_starts: str, exon_ends: str, cds_start: int, cds_end: int,
                 sv_start: int, sv_end: int):
    """completeCDSlength and overlappedCDSlength of one transcript."""
    pairs = list(itertools.zip_longest(exon_starts.split(","), exon_ends.split(","), fillvalue=""))
    complete = 0
    for a, b in pairs:
        if a == "" or b == "":
            continue
        a, b = int(a), int(b)
        if a < cds_start:
            if b < cds_start:
                continue
            complete += b - cds_start
            continue
        if b > cds_end:
            if a > cds_end:
                continue
            complete += cds_end - a
            continue
        complete += b - a
    overlapped = 0
    for a, b in pairs:
        if a == "" or b == "":
            continue
        a, b = int(a), int(b)
        if a < cds_start:
            if b < cds_start:
                continue
            a = cds_start
        if b > cds_end:
            if a > cds_end:
                continue
            b = cds_end
        if sv_start < a:
            if sv_end < a:
                break
            i = a
        else:
            i = sv_start
        if sv_start > b:
            continue
        j = sv_end if sv_end < b else b
        overlapped += j - i
    return complete, overlapped


def _selected_transcripts(config) -> List[str]:
    selected: List[str] = []
    if config["txFile"] != "":
        for line in lines_from_file(config["txFile"]):
            for tx in tcl.split_chars(line, " |\t"):
                selected.append(re.sub(r"\.[0-9]+", "", tx, count=1))
    return selected


def _format_input_bed(config, now: Optional[datetime.datetime] = None) -> str:
    """First part of genesAnnotation: trim, drop "chr", keep the standard
    contigs in karyotypic order, sort. Returns the new ``bedFile``."""
    c = config
    by_chrom: Dict[str, List[str]] = {}
    with open(c["bedFile"], "r", newline="") as f:
        for raw in f:
            line = raw[:-1] if raw.endswith("\n") else raw
            if line == "":
                continue
            line = re.sub("\t$", "\t.", line, count=1)
            line = "\t".join(v.strip() for v in line.split("\t"))
            first = tcl.lindex(line, 0)
            if "chr" in first:
                by_chrom.setdefault(first.replace("chr", "", 1), []).append(line[3:])
            else:
                by_chrom.setdefault(first, []).append(line)

    new_bed = re.sub(r".bed$", lambda _m: ".formatted.sorted.bed", c["bedFile"], count=1, flags=re.I)
    new_bed = os.path.join(c["outputDir"], os.path.basename(new_bed))
    if os.path.exists(new_bed):
        os.remove(new_bed)
    empty = True
    for chrom in KARYOTYPE:
        if chrom not in by_chrom:
            continue
        empty = False
        write_text_in_file("\n".join(by_chrom.pop(chrom)), f"{new_bed}.tmp")
    if empty:
        raise ConfigExit("The SV input file is empty or badly formatted.\nExit with error\n", 2)
    unix_sort_bed(f"{new_bed}.tmp", new_bed, "genesAnnotation")
    os.remove(f"{new_bed}.tmp")
    c["bedFile"] = new_bed
    return new_bed


def _bedtools_intersect(config, out_file: str) -> None:
    c = config
    base = [c.get("bedtools", "bedtools"), "intersect"]
    tail = ["-a", c["bedFile"], "-b", c["genesFile"], "-wa", "-wb"]
    for cmd in (base + ["-sorted"] + tail, base + tail):
        with open(out_file, "w") as out:
            r = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE)
        # Tcl "exec" fails on a non-zero exit status or on any stderr output
        if r.returncode == 0 and not r.stderr:
            return
    raise ConfigExit(f"-- genesAnnotation --\n{' '.join(base + ['-sorted'] + tail)} > {out_file}\n"
                     f"{r.stderr.decode(errors='replace')}\nExit with error", 2)


def genes_annotation(config) -> GenesResult:
    """genesAnnotation.

    Keeps 1 transcript per gene for each SV:
      - the one selected by the user with ``-txFile``
      - else the one with the most bp of overlapped CDS
      - else the one with the most bp of overlapped transcript.

    Writes ``<outputDir>/<outputFile>.tmp``: every SV line followed by "full",
    then one line per overlapped gene (SV line + gene file line +
    overlappedCDSlength, percentCDS, overlappedTxLength) followed by "split".
    """
    c = config
    check_svt_bed_col(c, c["bedFile"])
    check_samplesid_bed_col(c, c["bedFile"])
    _format_input_bed(c)
    bed_lines = lines_from_file(c["bedFile"])

    def bed_at(k: int) -> str:
        return bed_lines[k] if 0 <= k < len(bed_lines) else ""

    full_and_split = f"{c['outputDir']}/{c['outputFile']}.tmp"
    split_file = f"{full_and_split}.tmp"
    for p in (full_and_split, split_file):
        if os.path.exists(p):
            os.remove(p)
    result = GenesResult(full_and_split_file=full_and_split)

    _bedtools_intersect(c, split_file)
    if os.path.getsize(split_file) == 0:
        print("...no intersection between SV and gene annotation")
        write_text_in_file("\n".join(f"{l}\tfull" for l in bed_lines), full_and_split)
        os.remove(split_file)
        return result

    selected_tx = _selected_transcripts(c)
    to_write: List[str] = []
    all_genes: List[str] = []
    genes: List[str] = []
    annot: Dict[str, List[str]] = {}
    tx_len: Dict[str, List[int]] = {}
    cds_ov: Dict[str, List[int]] = {}
    cds_compl: Dict[str, List[int]] = {}
    finish: set = set()

    def write_split_lines(short_sv: str) -> None:
        nonlocal genes
        genes = remove_redundancy_without_sorting(genes)
        result.genes_by_sv[short_sv] = ";".join(genes)
        for gene in genes:
            best_ov, best_compl, best_tx, best_ann = -1, -1, -1, ""
            for ov, compl, txl, ann in zip(cds_ov[gene], cds_compl[gene], tx_len[gene], annot[gene]):
                if ov > best_ov or (ov == best_ov and txl > best_tx):
                    best_ann, best_ov, best_compl, best_tx = ann, ov, compl, txl
            percent = 0 if best_compl == 0 else best_ov * 100 // best_compl
            to_write.append(f"{best_ann}\t{best_ov}\t{percent}\t{best_tx}\tsplit")

    with open(split_file, "r", newline="") as f:
        split_lines = [l[:-1] if l.endswith("\n") else l for l in f]

    ls = split_lines[0].split("\t")
    sv_left, sv_right = int(ls[1]), int(ls[2])
    short_old = "\t".join(ls[0:3])
    old_split_sv = "\t".join(ls[:-10])

    n = 0
    sv_from_bed = bed_at(n)
    # Only 1 SV is associated with split lines AND it is not the first SV of the BED
    while sv_from_bed != old_split_sv:
        if n >= len(bed_lines):
            raise RuntimeError(f"genesAnnotation: SV not found in {c['bedFile']}: {old_split_sv}")
        to_write.append(f"{bed_at(n)}\tfull")
        n += 1
        sv_from_bed = bed_at(n)

    for line in split_lines:
        if line == "":
            continue
        ls = line.split("\t")
        split_sv = "\t".join(ls[:-10])
        if split_sv != old_split_sv:
            # all the split lines of the old SV are known: write them
            while sv_from_bed != split_sv:
                if n >= len(bed_lines):
                    raise RuntimeError(f"genesAnnotation: SV not found in {c['bedFile']}: {split_sv}")
                to_write.append(f"{bed_at(n)}\tfull")
                if sv_from_bed == old_split_sv:
                    write_split_lines(short_old)
                    sv_left, sv_right = int(ls[1]), int(ls[2])
                    annot, tx_len, cds_ov, cds_compl = {}, {}, {}, {}
                    finish = set()
                    genes = []
                    old_split_sv = split_sv
                    m = re.match("[^\t]+\t[^\t]+\t[^\t]+", old_split_sv)
                    if m:
                        short_old = m.group(0)
                n += 1
                sv_from_bed = bed_at(n)

        tx_name = ls[-5]
        gene = ls[-6]
        all_genes.append(gene)
        if gene in finish:
            continue
        genes.append(gene)
        annot.setdefault(gene, []).append(line)

        tx_start, tx_end = int(ls[-9]), int(ls[-8])
        cds_start, cds_end = int(ls[-4]), int(ls[-3])
        start = tx_start if tx_start > sv_left else sv_left
        end = tx_end if tx_end < sv_right else sv_right
        overlapped_tx = end - start
        tx_len.setdefault(gene, []).append(overlapped_tx)
        complete, overlapped_cds = _cds_lengths(ls[-2], ls[-1], cds_start, cds_end, start, end)
        cds_compl.setdefault(gene, []).append(complete)
        cds_ov.setdefault(gene, []).append(overlapped_cds)

        if tx_name != "" and tx_name in selected_tx:
            annot[gene] = [line]
            tx_len[gene] = [overlapped_tx]
            cds_ov[gene] = [overlapped_cds]
            cds_compl[gene] = [complete]
            finish.add(gene)

    # Last SV of the file
    to_write.append(f"{bed_at(n)}\tfull")
    n += 1
    write_split_lines(short_old)
    while n < len(bed_lines):
        to_write.append(f"{bed_at(n)}\tfull")
        n += 1

    write_text_in_file("\n".join(to_write), full_and_split)
    os.remove(split_file)
    result.all_genes_overlapped = sorted(set(all_genes))
    return result
