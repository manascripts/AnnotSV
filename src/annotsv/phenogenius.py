"""Port of AnnotSV-phenogenius.tcl (PhenoGenius, phenotype-driven gene ranking).

* ``check_phenogenius_cli``  -> ``checkPhenoGeniusCli``
* ``PhenoGenius.run``        -> ``runPhenoGeniusCli`` (writes and runs the same bash
  script: ``poetry install`` then ``poetry run python3 phenogenius_cli.py``)
* ``PhenoGenius.annotation`` -> ``PhenoGeniusCliAnnotation``

``config["PhenoGeniusCli"]`` is "1" or "0" like ``g_AnnotSV(PhenoGeniusCli)``.
"""

from __future__ import annotations

import os
import re
import time
from typing import Dict, List

from . import helpers, tcl
from .general import lines_from_file, write_text_in_file
from .ncbi import search_ncbi_gene_id


def _stamp() -> str:
    return f"{time.strftime('%Y%m%d-%H%M%S')}_{os.getpid()}"


def _rm(path: str) -> None:
    if os.path.isfile(path):
        os.remove(path)


def check_phenogenius_cli(config) -> None:
    """checkPhenoGeniusCli: PhenoGenius is used only with -hpo and a working install."""
    c = config
    if c["hpo"] == "":
        c["PhenoGeniusCli"] = "0"  # No HPO terms in INPUT
        return
    failed, message = helpers.check_phenogenius_cli_install(f"{c['installDir']}/share")
    if not failed:
        print("\tINFO: AnnotSV takes use of PhenoGenius (Yauy et al., 2023) for the phenotype-driven analysis.")
        c["PhenoGeniusCli"] = "1"
    else:
        log_file = f"{c['installDir']}/share/python3/phenogeniuscli/{_stamp()}.log"
        print("\nWARNING: No PhenoGenius installation available for the phenotype-driven analysis.")
        try:
            write_text_in_file(message, log_file)
        except OSError as exc:  # Tcl: error when the directory does not exist
            raise RuntimeError(f'couldn\'t open "{log_file}": {exc.strerror}') from exc
        print(f"\ncf {log_file}\n")
        c["PhenoGeniusCli"] = "0"


class PhenoGenius:
    """The ``g_PhenoGeniusCli`` array: gene name or NCBI ID -> "score\\tphenotype\\tspecificity"."""

    def __init__(self, config):
        self.config = config
        self.table: Dict[str, str] = {}

    @property
    def enabled(self) -> bool:
        return str(self.config.get("PhenoGeniusCli", "0")) == "1"

    def run(self, genes: List[str], ncbi_ids: List[str] = None) -> None:
        """runPhenoGeniusCli. ``ncbi_ids`` defaults to the NCBI IDs of ``genes``
        (as computed in regulatoryElementsAnnotation)."""
        c = self.config
        if ncbi_ids is None:
            ncbi_ids = sorted({i for i in (search_ncbi_gene_id(c, g) for g in genes) if i != ""})
        print(f"...running PhenoGeniusCli for {len(genes)} gene names "
              f"({time.strftime('%B %d %Y - %H:%M')})\n")
        stamp = _stamp()
        command_file = f"{c['outputDir']}/{stamp}_run_phenogeniuscli.tmp.bash"
        results_file = f"{c['outputDir']}/{stamp}_phenogeniuscli_results.tmp.tsv"
        log_file = f"{c['outputDir']}/{stamp}_phenogeniuscli_results.tmp.log"

        code = "#!/bin/bash\n\n"
        code += f"cd {c['installDir']}/share/python3/phenogeniuscli/PhenoGeniusCli\n"
        code += "poetry install\n"
        code += f"poetry run python3 phenogenius_cli.py --hpo_list {c['hpo']} --result_file {results_file} &> {log_file} \n"
        write_text_in_file(code, command_file)
        os.chmod(command_file, 0o755)

        failed, message = tcl.exec_command([command_file])
        if failed and re.search("error", message, re.I):
            print(f"{code}\n")
            print(message)

        if not os.path.exists(log_file):
            c["PhenoGeniusCli"] = "0"
            return

        self.parse(lines_from_file(results_file), genes, ncbi_ids)

        _rm(command_file)
        _rm(results_file)
        _rm(log_file)

    def parse(self, lines: List[str], genes: List[str], ncbi_ids: List[str]) -> None:
        genes_set, ids_set = set(genes), set(ncbi_ids)
        idx = None
        for L in lines:
            Ls = L.split("\t")
            # Header: #gene_id gene_symbol rank score hpo_implicated hpo_description_implicated phenotype_specificity
            if L.startswith("#gene_id"):
                idx = (0, _index(Ls, "gene_symbol"), _index(Ls, "score"),
                       _index(Ls, "hpo_description_implicated"), _index(Ls, "phenotype_specificity"))
                continue
            if idx is None:
                raise RuntimeError("can't read \"i_ncbi_id\": no such variable")
            i_id, i_gene, i_score, i_pheno, i_spec = idx
            ncbi_id = _li(Ls, i_id)
            gene = _li(Ls, i_gene)
            score = _li(Ls, i_score)
            phenotype = re.sub(r"': [0-9]\.[0-9]", "", _li(Ls, i_pheno))
            phenotype = re.sub(r"[\[\]{}']", "", phenotype)
            phenotype = phenotype.replace('"', "")
            specificity = _li(Ls, i_spec)[:1]  # A, B, C or D
            if gene in genes_set:
                self.table[gene] = f"{score}\t{phenotype}\t{specificity}"
            if ncbi_id in ids_set:
                self.table[ncbi_id] = f"{score}\t{phenotype}\t{specificity}"

    def annotation(self, gene: str, what: str) -> str:
        """PhenoGeniusCliAnnotation: what = "all" or "specificity"."""
        if what == "all":
            if gene in self.table:
                return self.table[gene]
            ncbi = search_ncbi_gene_id(self.config, gene)
            return self.table.get(ncbi, "-1.0\t\t")
        if what == "specificity":
            if gene in self.table:
                return self.table[gene].split("\t")[-1]
            ncbi = search_ncbi_gene_id(self.config, gene)
            if ncbi in self.table:
                return self.table[ncbi].split("\t")[-1]
            return ""
        print(f'proc PhenoGeniusCliAnnotation: Bad option value for "what" ({what})')
        return ""


def _li(ls: List[str], i: int) -> str:
    return ls[i] if 0 <= i < len(ls) else ""


def _index(ls: List[str], name: str) -> int:
    return ls.index(name) if name in ls else -1
