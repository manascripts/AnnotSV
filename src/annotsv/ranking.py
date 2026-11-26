"""SV ranking (ACMG classification): port of ``AnnotSV-ranking.tcl`` and of
the ranking part of ``OrganizeAnnotation`` (``AnnotSV-write.tcl``).

* ``Ranking.prepare``: SVprepareRanking (column positions, ``g_i``).
* ``Ranking.rank``: SVrankingLoss / SVrankingGain / SVrankingINS / SVrankingINV,
  called on every full and split line before it is written.
* ``Ranking.achieve``: achieveSVrankingLoss / achieveSVrankingGain.
* ``Ranking.sort_ids``: the output order (ranking score, then Exomiser score,
  then input order).
* ``Ranking.full_columns``: AnnotSV_ranking_score, AnnotSV_ranking_criteria
  and ACMG_class of a full line.

The Exomiser and PhenoGenius scores are read from the
``Exomiser_gene_pheno_score`` / ``PhenoGenius_specificity`` columns when the
header has them. The best Exomiser score per SV (the second sort key) comes
from the Exomiser module, not ported yet: ``sort_ids`` takes it as an argument
and defaults to ``-1`` like Tcl does when Exomiser is off.

Values are compared with Tcl ``expr`` semantics (``tcl.compare``): numeric
when both sides look like numbers, string comparison otherwise. So an empty
``DDD_HI_percent`` counts as ``<= 10`` and an empty ``Overlapped_CDS_percent``
as ``< 10`` (Tcl parity).
"""

from __future__ import annotations

import functools
import math
import re
from typing import Dict, Iterable, List, Optional, Tuple

from . import tcl
from .config import DEFAULT_RANKING_CRITERIA
from .general import normalize_svtype

# SVprepareRanking: (key, header regexp, compulsory)
_COLUMNS = [
    ("gene", "Gene_name", True), ("NbGenes", "Gene_count", True), ("full", "Annotation_mode", True),
    ("REgene", "RE_gene", True),
    ("Ploss", "P_loss_coord", True), ("Pgain", "P_gain_coord", True), ("Bloss", "B_loss_coord", True),
    ("Bgain", "B_gain_coord", True), ("Psnvindel", "P_snvindel_nb", True),
    ("poPloss", "po_P_loss_coord", True), ("poPgain", "po_P_gain_coord", True),
    ("poBlossAllG", "po_B_loss_allG_coord", True), ("poBlossSomeG", "po_B_loss_someG_coord", True),
    ("poBgainAllG", "po_B_gain_allG_coord", True), ("poBgainSomeG", "po_B_gain_someG_coord", True),
    ("HI", "HI", True), ("TS", "TS", True), ("pLI", "GnomAD_pLI", True), ("loeuf", "LOEUF_bin", True),
    ("HIpercent", "DDD_HI_percent", True),
    ("exomiser", "Exomiser_gene_pheno_score", False), ("phenogenius", "PhenoGenius_specificity", False),
    ("morbid", "OMIM_morbid", True),
    ("location", "Location", True), ("location2", "Location2", True),
    ("CDSlength", "Overlapped_CDS_length", True), ("CDSpercent", "Overlapped_CDS_percent", True),
    ("frameshift", "Frameshift", True), ("NbExons", "Exon_count", True),
]

# achieveSVrankingLoss / achieveSVrankingGain: section 2, first match only (highest score first)
_LOSS_S2 = ["2A", "2C-1", "2D-2", "2D-4", "2E-1", "2E-2", "2D-3", "2E-3", "2E-4", "2H", "2F"]
_GAIN_S2 = ["2A", "2I-1", "2I-2", "2K", "2H-1", "2C", "2D", "2F"]
_GENE_LIST_LOSS = {"2C-1", "2D-2", "2D-4", "2E-1", "2E-2", "2D-3", "2E-3", "2H"}
_GENE_LIST_GAIN = {"2I-1", "2I-2", "2K", "2H-1"}
_RE_EX = re.compile(r"(.+?) \(.*?EX=(.+)$", re.S)


def _ge(a, b) -> bool:
    return tcl.compare(str(a), str(b)) >= 0


def _gt(a, b) -> bool:
    return tcl.compare(str(a), str(b)) > 0


def _le(a, b) -> bool:
    return tcl.compare(str(a), str(b)) <= 0


def _lt(a, b) -> bool:
    return tcl.compare(str(a), str(b)) < 0


def _num(value):
    n = tcl.to_number(value) if isinstance(value, str) else value
    if n is None:
        raise ValueError(f'can\'t use non-numeric string "{value}" as operand of "+"')
    return n


def _add(a, b):
    """``expr {$a+$b}``: integer when both are integers, double otherwise."""
    return _num(a) + _num(b)


def tcl_str(value) -> str:
    """How Tcl prints an ``expr`` result (shortest round-trip double, like ``repr``)."""
    if isinstance(value, float):
        if math.isinf(value):
            return "Inf" if value > 0 else "-Inf"
        r = repr(value)
        return r
    return str(value)


def _tcl_round(x) -> int:
    """Tcl ``round()``: half away from zero."""
    if isinstance(x, int):
        return x
    if x >= 0:
        f = math.floor(x)
        return int(f) + (1 if x - f >= 0.5 else 0)
    f = math.floor(-x)
    return -(int(f) + (1 if -x - f >= 0.5 else 0))


def round_score(score) -> float:
    """``double(round(100*$score))/100`` (e.g. 2.1999999999999997 -> 2.2)."""
    return _tcl_round(100 * _num(score)) / 100


def acmg_class(score) -> str:
    if score == "":
        return "NA"
    if _ge(tcl_str(score), "0.99"):
        return "5"
    if _ge(tcl_str(score), "0.9"):
        return "4"
    if _ge(tcl_str(score), "-0.89"):
        return "3"
    if _ge(tcl_str(score), "-0.98"):
        return "2"
    return "1"


def _join_genes(genes: List[str]) -> str:
    return "/".join(sorted(set(genes)))


def _best_re_gene(re_gene: str) -> Tuple[str, Optional[str]]:
    """Gene of RE_gene with the best "EX=" (Exomiser) score."""
    best_ex, best_g = "0", None
    for infos in tcl.split_chars(re_gene, "\\);"):
        m = _RE_EX.search(infos)
        if m:
            g, ex = m.group(1), m.group(2).replace(",", ".", 1)
            if _gt(ex, best_ex):
                best_ex, best_g = ex, g
    return best_ex, best_g


class Ranking:
    """The ``g_i``, ``g_rankingScore`` and ``g_rankingExplanations`` globals."""

    def __init__(self, config):
        self.c = config
        crit = dict(DEFAULT_RANKING_CRITERIA)
        crit.update(getattr(config, "ranking_criteria", {}) or {})
        self.crit = crit
        self.i: Dict[str, int] = {}
        self.enabled = False
        self.score: Dict[str, object] = {}
        self.expl: Dict[str, Dict[str, object]] = {}   # AnnotSV_ID -> {"": text, "2A": text, "5H": [genes]}
        self.full_done: set = set()
        self.split_done: set = set()
        self.achieved: set = set()
        self.acmg_full: Dict[str, str] = {}

    # ------------------------------------------------------------------
    # SVprepareRanking
    # ------------------------------------------------------------------
    def prepare(self, header: str) -> bool:
        c = self.c
        self.enabled = False
        c["ranking"] = "0"
        ls = header.split("\t")
        # VCF input: SV_type is the 6th column of the output
        if re.search(".vcf$|.vcf.gz$", str(c["SVinputFile"]), re.I):
            c["svtTSVcol"] = "5"
            c["svtBEDcol"] = "3"
        if str(c.get("svtTSVcol", "-1")) == "-1":
            self._warning(header, "SV_type")
            return False
        self.i = {}
        for key, pattern, compulsory in _COLUMNS:
            idx = next((k for k, col in enumerate(ls) if re.search(pattern, col)), -1)
            if idx == -1 and compulsory:
                self.i = {}
                self._warning(header, pattern)
                return False
            self.i[key] = idx
        self.enabled = True
        c["ranking"] = "1"
        return True

    def _warning(self, header: str, field: str) -> None:
        if self.c.get("organism") == "Human":
            print(f'WARNING: "{field}" annotation is missing:')
            print("******************************************")
            print(header)
            print("******************************************")
            print("         => No ranking provided\n")

    # ------------------------------------------------------------------
    # SVranking*
    # ------------------------------------------------------------------
    def rank(self, text: str, svtype: str) -> None:
        """Called on every output line (full and split), ``svtype`` from the input."""
        svtype = normalize_svtype(svtype)
        if svtype == "DEL":
            self._loss(text)
        elif svtype == "DUP":
            self._gain(text)
        else:   # INS, INV and the others: no ranking
            aid = text.split("\t")[0]
            self.score[aid] = ""
            self.expl[aid] = {"": ""}

    def _col(self, ls: List[str], key: str) -> str:
        k = self.i[key]
        return ls[k] if 0 <= k < len(ls) else ""

    def _e(self, aid: str) -> Dict[str, object]:
        return self.expl.setdefault(aid, {})

    def _append(self, aid: str, key: str, gene: str) -> None:
        self._e(aid).setdefault(key, []).append(gene)

    def _section1(self, aid, re_gene, nb_genes, kind):
        cr = self.crit
        e = self._e(aid)
        if re_gene == "" and nb_genes == "0":
            self.score[aid] = cr[f"1B_{kind}"]
            e[""] = f"1B (cf Gene_count, RE_gene, {cr[f'1B_{kind}']});"
        else:
            self.score[aid] = cr[f"1A_{kind}"]
            e[""] = f"1A (cf Gene_count, RE_gene, +{cr[f'1A_{kind}']});"

    def _section3(self, aid, nb_genes, kind, max_a, max_b):
        cr = self.crit
        e = self._e(aid)
        if _le(nb_genes, max_a):
            word = "gene" if _lt(nb_genes, 2) else "genes"
            e["3A"] = f"3A ({nb_genes} {word}, +{cr[f'3A_{kind}']});"
        elif _le(nb_genes, max_b):
            e["3B"] = f"3B ({nb_genes} genes, +{cr[f'3B_{kind}']});"
        else:
            e["3C"] = f"3C ({nb_genes} genes, +{cr[f'3C_{kind}']});"

    def _section5_re(self, aid, re_gene):
        if re_gene != "":
            best_ex, best_g = _best_re_gene(re_gene)
            if _ge(best_ex, "0.7"):
                self._append(aid, "5H", f"RE:{best_g}")
            elif _ge(best_ex, "0.5"):
                self._append(aid, "5G", f"RE:{best_g}")

    def _section5_split(self, aid, exomiser, phenogenius, gene):
        if _ge(exomiser, "0.7") or phenogenius == "A":
            self._append(aid, "5H", gene)
        elif _ge(exomiser, "0.5") or phenogenius == "B":
            self._append(aid, "5G", gene)

    def _loss(self, text: str) -> None:
        cr = self.crit
        ls = text.split("\t")
        aid = ls[0]
        col = lambda k: self._col(ls, k)  # noqa: E731
        if col("full") == "full":
            # In case of SV redundancy in the input file
            if aid in self.full_done:
                self.split_done.add(aid)
                return
            self.full_done.add(aid)
            nb_genes, re_gene = col("NbGenes"), col("REgene")
            self._section1(aid, re_gene, nb_genes, "Loss")
            e = self._e(aid)
            # Section 2
            if col("Ploss") != "":
                e["2A"] = f"2A (cf P_loss_source, +{cr['2A_Loss']});"
            else:
                if col("poPloss") != "":
                    e["2B"] = f"2B (cf po_P_loss_source, HI, OMIM_morbid, +{cr['2B_Loss']});"
                if col("Bloss") != "":
                    e["2F"] = f"2F/4O (cf B_loss_source, {cr['2F-4O_Loss']});"
                elif col("poBlossSomeG") != "":
                    e["2G"] = f"2G (cf po_B_loss_someG_source, +{cr['2G_Loss']});"
            # Section 3
            self._section3(aid, nb_genes, "Loss", "24", "35")
            # Section 4
            if col("poBlossAllG") != "" and "2F" not in e:
                e["4O"] = f"4O (cf B_loss_source, po_B_loss_allG_source, {cr['4O_Loss']});"
            # Section 5
            self._section5_re(aid, re_gene)
            return

        # Split lines
        if aid in self.split_done:
            return
        hi, morbid, gene = col("HI"), col("morbid"), col("gene")
        exomiser = col("exomiser").replace(",", ".", 1)
        phenogenius = col("phenogenius")
        if hi == "3" or morbid == "yes":
            location, location2 = col("location"), col("location2")
            cds_length = col("CDSlength")
            cds_percent = col("CDSpercent").replace(",", ".", 1)
            psnvindel, nb_exons = col("Psnvindel"), col("NbExons")
            if re.search("^5'UTR", location2) and not re.search("3'UTR$", location2):
                self._append(aid, "2C-1" if cds_length != "0" else "2C-2", gene)
            elif re.search("3'UTR$", location2) and not re.search("^5'UTR", location2):
                last_exon = f"exon{nb_exons}-txEnd"
                if cds_length == "0":
                    self._append(aid, "2D-1", gene)
                elif location == last_exon and psnvindel != "":
                    self._append(aid, "2D-2", gene)
                elif location == last_exon:
                    self._append(aid, "2D-3", gene)
                elif location != last_exon and location != f"intron{tcl_str(_add(nb_exons, -1))}-txEnd":
                    self._append(aid, "2D-4", gene)
            elif not re.search("txStart|txEnd", location):
                frameshift = col("frameshift")
                exon_overlapped = False
                if "exon" in location:
                    exon_overlapped = True
                elif "intron" in location:
                    # Tcl parity: compares the characters 0 and 2 of "i-j" (single digit introns)
                    tutu = location.replace("intron", "")
                    if tutu[0:1] != tutu[2:3]:
                        exon_overlapped = True
                if frameshift == "yes":
                    self._append(aid, "2E-1", gene)
                elif exon_overlapped and psnvindel != "" and _ge(cds_percent, "10"):
                    self._append(aid, "2E-2", gene)
                elif exon_overlapped and psnvindel != "" and _lt(cds_percent, "10"):
                    self._append(aid, "2E-3", gene)
                elif exon_overlapped and psnvindel == "" and _gt(cds_percent, "10"):
                    self._append(aid, "2E-4", gene)
        # 2H
        pli = col("pLI").replace(",", ".")
        hi_percent = col("HIpercent").replace(",", ".")
        n = (1 if _ge(pli, "0.9") else 0) + (1 if _le(hi_percent, "10") else 0)
        if n >= 2:
            self._append(aid, "2H", gene)
        # Section 5
        self._section5_split(aid, exomiser, phenogenius, gene)

    def _gain(self, text: str) -> None:
        cr = self.crit
        ls = text.split("\t")
        aid = ls[0]
        col = lambda k: self._col(ls, k)  # noqa: E731
        if col("full") == "full":
            nb_genes, re_gene = col("NbGenes"), col("REgene")
            pgain, bgain = col("Pgain"), col("Bgain")
            po_pgain, po_some, po_all = col("poPgain"), col("poBgainSomeG"), col("poBgainAllG")
            self._section1(aid, re_gene, nb_genes, "Gain")
            e = self._e(aid)
            # Section 2 (Tcl parity: 2A and 2B write "+1.00" / "+0.00", not the criteria)
            if pgain != "":
                e["2A"] = "2A (cf P_gain_source, +1.00);"
            else:
                if po_pgain != "":
                    e["2B"] = "2B (cf po_P_gain_source, TS, OMIM_morbid, +0.00);"
                if po_all != "":
                    e["2C"] = f"2C (cf po_B_gain_allG_source, {cr['2C_Gain']});"
                elif bgain != "":
                    e["2DE"] = ""
                elif po_some != "":
                    e["2G"] = f"2G (cf po_B_gain_someG_source, +{cr['2G_Gain']});"
            # Section 3
            self._section3(aid, nb_genes, "Gain", "34", "49")
            # Section 4
            if po_all != "" or bgain != "" or po_some != "":
                e["4O"] = (f"4O (cf po_B_gain_AllG_source, B_gain_source, po_B_gain_SomeG_source, "
                           f"+{cr['4O_Gain']});")
            # Section 5
            self._section5_re(aid, re_gene)
            return

        # Split lines
        hi, morbid, gene, location = col("HI"), col("morbid"), col("gene"), col("location")
        exomiser = col("exomiser").replace(",", ".", 1)
        phenogenius = col("phenogenius")
        e = self._e(aid)
        if "2DE" in e:
            if not re.search("exon|intron", location):
                e["2D"] = f"2D (cf B_gain_source, {cr['2D_Gain']});"
            else:
                e["2E"] = f"2E (cf B_gain_source, +{cr['2E_Gain']});"
        pheno_a = _ge(exomiser, "0.7") or phenogenius == "A"
        starts, ends = location.startswith("txStart"), location.endswith("txEnd")
        if hi == "3" or morbid == "yes":
            if not starts and not ends:
                if col("frameshift") == "yes":
                    self._append(aid, "2I-1", gene)
                elif pheno_a:
                    self._append(aid, "2I-2", gene)
                else:
                    self._append(aid, "2I-3", gene)
            elif not starts or not ends:
                self._append(aid, "2K" if pheno_a else "2J", gene)
            elif location == "txStart-txEnd":
                self._append(aid, "2H-1" if pheno_a else "2H-2", gene)
        elif not starts or not ends:
            self._append(aid, "2L", gene)
        self._section5_split(aid, exomiser, phenogenius, gene)

    # ------------------------------------------------------------------
    # achieveSVranking*
    # ------------------------------------------------------------------
    def achieve_all(self, ids: Iterable[str]) -> None:
        for aid in sorted(set(ids)):
            svtype = normalize_svtype(aid.split("_")[-2] if len(aid.split("_")) >= 2 else aid)
            if svtype == "DEL":
                self._achieve_loss(aid)
            elif svtype == "DUP":
                self._achieve_gain(aid)

    def _add_score(self, aid, value):
        self.score[aid] = _add(self.score[aid], value)

    def _genes_text(self, aid, key, crit_key, sep=", +", end=");"):
        return f"{key} ({_join_genes(self.expl[aid][key])}{sep}{self.crit[crit_key]}{end}"

    def _achieve_loss(self, aid: str) -> None:
        cr = self.crit
        if aid in self.achieved:
            return
        self.achieved.add(aid)
        if aid not in self.score:
            self.score[aid] = ""
            return
        e = self._e(aid)
        e.setdefault("", "")
        # Section 2
        for key in _LOSS_S2:
            if key in e:
                ck = "2F-4O_Loss" if key == "2F" else f"{key}_Loss"
                self._add_score(aid, cr[ck])
                if key == "2E-4":
                    # Tcl parity: "2E-4 (genes), +0.20;"
                    e[""] += self._genes_text(aid, key, ck, sep="), +", end=";")
                elif key in _GENE_LIST_LOSS:
                    e[""] += self._genes_text(aid, key, ck)
                else:
                    e[""] += e[key]
                break
        if "2B" in e and not re.search("2A|2C|2D|2E", e[""]):
            e[""] += e["2B"]
        if "2C-2" in e:
            e[""] += self._genes_text(aid, "2C-2", "2C-2_Loss", sep=",+")
        if "2D-1" in e:
            e[""] += self._genes_text(aid, "2D-1", "2D-1_Loss", sep=",+")
        if "2G" in e:
            e[""] += e["2G"]
        # Section 3
        for key in ("3C", "3B", "3A"):
            if key in e:
                self._add_score(aid, cr[f"{key}_Loss"])
                e[""] += e[key]
                break
        # Section 4
        if "4O" in e and not any(k in e for k in ("2A", "2F", "2C-1", "2D-2", "2D-3", "2D-4", "2E-1", "2E-2",
                                                  "2E-3", "2E-4", "2EH")):
            self._add_score(aid, cr["4O_Loss"])
            e[""] += e["4O"]
        # Section 5 (Tcl parity: no "+" before the 5H criteria)
        if "5H" in e:
            self._add_score(aid, cr["5H_Loss"])
            e[""] += f"5H ({_join_genes(e['5H'])}, {cr['5H_Loss']})"
        elif "5G" in e:
            self._add_score(aid, cr["5G_Loss"])
            e[""] += f"5G ({_join_genes(e['5G'])}, +{cr['5G_Loss']})"
        else:
            self._add_score(aid, cr["5F_Loss"])
            e[""] += f"5F (+{cr['5F_Loss']})"

    def _achieve_gain(self, aid: str) -> None:
        cr = self.crit
        if aid not in self.score:
            self.score[aid] = ""
            return
        e = self._e(aid)
        e.setdefault("", "")
        # Section 2: first match only
        for key in _GAIN_S2:
            if key in e:
                ck = f"{key}_Gain"
                self._add_score(aid, cr[ck])
                if key in _GENE_LIST_GAIN:
                    e[""] += self._genes_text(aid, key, ck)
                else:
                    e[""] += e[key]
                break
        # The "0" scores: all reported
        for key in ("2B", "2E", "2G"):
            if key in e:
                self._add_score(aid, cr[f"{key}_Gain"])
                e[""] += e[key]
        for key in ("2H-2", "2I-3", "2J", "2L"):
            if key in e:
                self._add_score(aid, cr[f"{key}_Gain"])
                e[""] += self._genes_text(aid, key, f"{key}_Gain")
        # Section 3
        for key in ("3C", "3B", "3A"):
            if key in e:
                self._add_score(aid, cr[f"{key}_Gain"])
                e[""] += e[key]
                break
        # Section 4
        if "4O" in e:
            self._add_score(aid, cr["4O_Gain"])
            e[""] += e["4O"]
        # Section 5
        if "2K" not in e:
            if "5H" in e:
                # Tcl parity: adds 0.15 whatever the 5H_Gain criteria
                self._add_score(aid, "0.15")
                e[""] += f"5H ({_join_genes(e['5H'])}, +{cr['5H_Gain']})"
            elif "5G" in e:
                self._add_score(aid, cr["5G_Gain"])
                e[""] += f"5G ({_join_genes(e['5G'])}, +{cr['5G_Gain']})"
            else:
                self._add_score(aid, cr["5F_Gain"])
                e[""] += f"5F (+{cr['5F_Gain']})"

    # ------------------------------------------------------------------
    # OrganizeAnnotation: order and ranking columns
    # ------------------------------------------------------------------
    def sort_ids(self, ids: List[str], exomiser_scores: Optional[Dict[str, str]] = None) -> List[str]:
        """Higher ranking score first, then higher Exomiser score, then input order."""
        exomiser_scores = exomiser_scores or {}
        rows = []
        for aid in ids:
            s = self.score.get(aid, "")
            rows.append((aid, "-99" if s == "" else tcl_str(s), exomiser_scores.get(aid, "-1")))
        # lsort -command DescendingSortOnElement1 [lsort -command DescendingSortOnElement2 ...]:
        # stable merge sorts with expr comparisons
        for k in (2, 1):
            rows = sorted(rows, key=functools.cmp_to_key(lambda x, y, k=k: -tcl.compare(x[k], y[k])))
        return [r[0] for r in rows]

    def full_columns(self, aid: str, metrics_fr: bool) -> Tuple[str, str, str]:
        """Ranking score, criteria and ACMG class of a full line (also kept for the split lines)."""
        score = self.score.get(aid, "")
        if score != "":
            score = round_score(score)
            self.score[aid] = score
        score_text = tcl_str(score)
        if metrics_fr:
            score_text = score_text.replace(".", ",")
        e = self.expl.get(aid)
        if e is None or "" not in e:
            criteria = ""
        else:
            if metrics_fr:
                e[""] = re.sub(r"([0-9])\.([0-9])", r"\1,\2", e[""])
            criteria = e[""]
        cls = acmg_class(score)
        self.acmg_full[aid] = f"full={cls}"
        return score_text, criteria, cls
