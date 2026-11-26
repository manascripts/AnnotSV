"""Port of ``searchforNCBIGeneID`` (AnnotSV-ncbi.tcl).

The ``geneID`` Tcl array is read once from
``Gene-based/NCBIandHGNCgeneID/results.txt`` and kept on the config object.

Tcl quirk kept: the header columns are searched as ``"Previous symbol"`` and
``"Alias symbol"`` (singular), while the HGNC file says ``"Previous symbols"``
and ``"Alias symbols"``. ``lsearch`` then gives -1, ``lindex`` gives "" and
every line sets ``geneID("")``: the empty gene name gets the NCBI ID of the
last line.
"""

from __future__ import annotations

import re
from typing import Dict, List

from .general import lines_from_file


def _li(ls: List[str], i: int) -> str:
    return ls[i] if 0 <= i < len(ls) else ""


def _index(ls: List[str], name: str) -> int:
    return ls.index(name) if name in ls else -1


def _load(config) -> Dict[str, str]:
    c = config
    path = f"{c['annotationsDir']}/Annotations_{c['organism']}/Gene-based/NCBIandHGNCgeneID/results.txt"
    gene_id: Dict[str, str] = {}
    idx = None
    for L in lines_from_file(path):
        Ls = L.split("\t")
        if re.search("Approved symbol", L):
            idx = (_index(Ls, "Approved symbol"), _index(Ls, "NCBI Gene ID"),
                   _index(Ls, "Previous symbol"), _index(Ls, "Alias symbol"))
            continue
        if idx is None:
            raise RuntimeError("can't read \"i_approuved\": no such variable")
        i_app, i_ncbi, i_prev, i_alias = idx
        ncbi = _li(Ls, i_ncbi)
        if "[" in ncbi:  # some lines have a bad "NCBI_gene_ID"
            continue
        gene_id[_li(Ls, i_app)] = ncbi
        gene_id[_li(Ls, i_prev)] = ncbi
        gene_id[_li(Ls, i_alias)] = ncbi
    return gene_id


def search_ncbi_gene_id(config, gene_name: str) -> str:
    """searchforNCBIGeneID: the NCBI gene ID of a gene name (or "")."""
    table = getattr(config, "_ncbi_gene_id", None)
    if table is None:
        table = _load(config)
        config._ncbi_gene_id = table
    return table.setdefault(gene_name, "")
