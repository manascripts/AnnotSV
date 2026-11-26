"""Port of AnnotSV-exomiser.tcl (Exomiser, phenotype-driven gene ranking).

* ``check_exomiser_installation`` -> ``checkExomiserInstallation``
* ``start_rest_service``          -> ``startTheRESTservice`` (``java -jar`` in the background)
* ``Exomiser.run``                -> ``runExomiser`` (one HTTP request per gene)
* ``Exomiser.annotation``         -> ``ExomiserAnnotation``

``config["Exomiser"]`` is "1" or "0" like ``g_AnnotSV(Exomiser)``, and
``config["hpoVersion"]`` stands for the ``hpoVersion`` Tcl global.

The JSON numbers are kept as their literal text (Tcl ``json2dict``), JSON
null/true/false become "null"/"true"/"false". The Tcl HTTP client ignores
the proxy environment variables, so the requests here bypass them too.
"""

from __future__ import annotations

import glob
import json
import os
import re
import signal
import socket
import subprocess
import time
import urllib.error
import urllib.request
from typing import Dict, List, Optional

from . import tcl
from .general import content_from_file, lines_from_file, write_text_in_file
from .ncbi import search_ncbi_gene_id

# startTheRESTservice: "after 6000" between two connection attempts, 50 attempts max (5 min)
WAIT_SECONDS = 6.0
MAX_ATTEMPTS = 50


def _stamp() -> str:
    return f"{time.strftime('%Y%m%d-%H%M%S')}_{os.getpid()}"


def _rm(path: str) -> None:
    if os.path.isfile(path):
        os.remove(path)


def check_exomiser_installation(config) -> None:
    """checkExomiserInstallation: Exomiser is used only with -hpo and its data installed."""
    c = config
    if c["hpo"] != "":
        c["Exomiser"] = "1"
    else:
        c["Exomiser"] = "0"
        return
    exo_dir = f"{c['annotationsDir']}/Annotations_Exomiser"
    ok = [os.path.basename(p) for p in glob.glob(f"{exo_dir}/*")
          if re.search(r"^([0-9]+)$", os.path.basename(p))]
    if ok:
        hpo_version = sorted(ok, key=int)[-1]
        c["hpoVersion"] = hpo_version
        print("\tINFO: AnnotSV takes use of Exomiser (Smedley et al., 2015) for the phenotype-driven analysis.")
        print(f"\tINFO: AnnotSV is using the Human Phenotype Ontology (version {hpo_version}). "
              f"Find out more at http://www.human-phenotype-ontology.org")
        if not os.path.exists(f"{c['etcDir']}/application.properties") and \
                not os.path.exists(f"{exo_dir}/{hpo_version}/application.properties"):
            print("\nWARNING: No Exomiser application properties file available")
            print(f"...{c['etcDir']}/application.properties doesn't exist\n")
            print(f"...{exo_dir}/{hpo_version}/application.properties doesn't exist\n")
            c["Exomiser"] = "0"
    else:
        print(f"\nWARNING: No Exomiser annotations available in {exo_dir}/\n")
        c["Exomiser"] = "0"


def _connect(port: int) -> bool:
    try:
        with socket.create_connection(("localhost", port), timeout=10):
            return True
    except OSError:
        return False


def start_rest_service(config, properties_file: str, port: int, start_file: str) -> Optional[subprocess.Popen]:
    """startTheRESTservice: the java process (or None if the service has not been started)."""
    c = config
    jars = glob.glob(f"{c['annotationsDir']}/jar/exomiser-rest-prioritiser-*.jar")
    jar = jars[-1] if jars else ""
    cmd = ["java", "-Xmx4g", "-jar"] + ([jar] if jar else []) + \
          [f"--server.port={port}", f"--spring.config.location={properties_file}"]
    text_cmd = f"java -Xmx4g -jar {jar} --server.port={port} --spring.config.location={properties_file} >& {start_file}"
    try:
        with open(start_file, "wb") as out:
            proc = subprocess.Popen(cmd, stdin=subprocess.DEVNULL, stdout=out, stderr=subprocess.STDOUT)
    except OSError as exc:
        print("\nWARNING: No Exomiser annotations available.")
        print("The REST service has not been started successfully:")
        print(text_cmd)
        print(f'couldn\'t execute "java": {(exc.strerror or str(exc)).lower()}')
        c["Exomiser"] = "0"
        _rm(start_file)
        return None

    # Wait (max 5 minutes) for the start
    waiting, i = True, 1
    while waiting:
        time.sleep(WAIT_SECONDS)
        if not _connect(port):
            i += 1
            if i == MAX_ATTEMPTS:
                waiting = False
        else:
            waiting = False
        if os.path.exists(start_file):
            content = content_from_file(start_file)
            if re.search("APPLICATION FAILED TO START|Application run failed", content):
                i, waiting = MAX_ATTEMPTS, False
            for L in lines_from_file(start_file):
                if L.startswith("ERROR"):
                    print(L)
                    i, waiting = MAX_ATTEMPTS, False

    if i == MAX_ATTEMPTS:
        print("\nWARNING: No Exomiser annotations available.")
        print("The REST service has not been started successfully:")
        print(text_cmd)
        print(f"(see {start_file})")
        c["Exomiser"] = "0"
        # Tcl does not kill the java process here: it is left running
        return None
    print(f"\t...idService = {proc.pid}\n")
    return proc


def _tcl_value(v) -> str:
    """A JSON value as the text Tcl json2dict gives (numbers already kept as text)."""
    if v is None:
        return "null"
    if v is True:
        return "true"
    if v is False:
        return "false"
    if isinstance(v, str):
        return v
    if isinstance(v, list):
        return tcl.format_list([_tcl_value(x) for x in v])
    return tcl.format_list([x for k, val in v.items() for x in (k, _tcl_value(val))])


def _get(d, *keys):
    """``dict get $d k1 k2...``: KeyError/TypeError when it fails."""
    for k in keys:
        if not isinstance(d, dict):
            raise TypeError("missing value to go with key")
        d = d[k]
    return d


def parse_exomiser_result(text: str) -> Optional[str]:
    """One REST answer -> "score\\thuman\\tmouse\\tfish" (or None: the gene is skipped)."""
    d_all = json.loads(text, parse_float=str, parse_int=str)
    try:
        results = _get(d_all, "results")
        if not isinstance(results, list) or not results:
            raise TypeError
        d_results = results[0]
        score = _get(d_results, "score")
    except (KeyError, TypeError):
        return None
    score = f"{float(_tcl_value(score)):.4f}"
    try:
        evidence = _get(d_results, "phenotypeEvidence")
    except (KeyError, TypeError):
        return None
    pheno: Dict[str, List[str]] = {"HUMAN": [], "MOUSE": [], "FISH": []}
    for d_organism in (evidence if isinstance(evidence, list) else []):
        try:
            organism = _tcl_value(_get(d_organism, "model", "organism"))
        except (KeyError, TypeError):
            continue
        acc = pheno.setdefault(organism, [])
        try:
            acc.append(_tcl_value(_get(d_organism, "model", "diseaseTerm")))
        except (KeyError, TypeError):
            pass
        try:
            best = _get(d_organism, "bestModelPhenotypeMatches")
        except (KeyError, TypeError):
            continue
        for d_query in (best if isinstance(best, list) else []):
            for k in ("query", "match"):
                try:
                    acc.append(_tcl_value(_get(d_query, k, "label")))
                except (KeyError, TypeError):
                    pass
    return "\t".join([score] + [";".join(sorted(set(pheno[o]))) for o in ("HUMAN", "MOUSE", "FISH")])


class Exomiser:
    """The ``g_Exomiser`` array: gene name -> "score\\thuman\\tmouse\\tfish" evidence."""

    def __init__(self, config):
        self.config = config
        self.table: Dict[str, str] = {}

    @property
    def enabled(self) -> bool:
        return str(self.config.get("Exomiser", "0")) == "1"

    def _free_port(self) -> int:
        failed, message = tcl.exec_command(["bash", f"{self.config['bashDir']}/searchForAFreePortNumber.sh"])
        if not failed:
            try:
                return int(message)
            except ValueError:
                pass
        print(message)
        print("WARNING: port is defined to 50000")
        return 50000

    def run(self, genes: List[str]) -> None:
        """runExomiser."""
        c = self.config
        print(f"...running Exomiser on {len(genes)} gene names ({time.strftime('%B %d %Y - %H:%M')})")
        port = self._free_port()
        print(f"\t...on port {port}")

        hpo_version = c.get("hpoVersion", "")
        exo_dir = f"{c['annotationsDir']}/Annotations_Exomiser/{hpo_version}"
        properties_file = f"{c['outputDir']}/{_stamp()}_exomiser_application.properties"
        if os.path.exists(f"{c['etcDir']}/application.properties"):
            infos = content_from_file(f"{c['etcDir']}/application.properties")
        else:
            infos = content_from_file(f"{exo_dir}/application.properties")
        infos = re.sub("XXXX", lambda _m: str(port), infos, count=1)
        infos = re.sub("YYYY", lambda _m: exo_dir, infos, count=1)
        infos = re.sub("ZZZZ", lambda _m: hpo_version, infos, count=1)
        write_text_in_file(infos, properties_file)

        start_file = f"{c['outputDir']}/{_stamp()}_exomiser.tmp"
        print("\t...starting the REST service")
        proc = start_rest_service(c, properties_file, port, start_file)
        if proc is not None:
            opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
            for gene in genes:
                gene_id = search_ncbi_gene_id(c, gene)
                if gene_id == "":
                    continue
                url = (f"http://localhost:{port}/exomiser/api/prioritise/?phenotypes={c['hpo']}"
                       f"&prioritiser=hiphive&prioritiser-params=human,mouse,fish,ppi&genes={gene_id}")
                try:
                    try:
                        with opener.open(url) as r:
                            body, charset = r.read(), r.headers.get_content_charset()
                    except urllib.error.HTTPError as e:  # Tcl reads the body whatever the status
                        body, charset = e.read(), e.headers.get_content_charset()
                    text = body.decode(charset or "utf-8", errors="replace")
                    result = parse_exomiser_result(text)
                except Exception as exc:  # noqa: BLE001 (Tcl "catch")
                    print(f"geneName: {gene}")
                    print(f"{exc}\n")
                    print(url)
                    continue
                if result is not None:
                    self.table[gene] = result
            # End the REST service
            try:
                os.kill(proc.pid, signal.SIGKILL)
                proc.wait()
            except OSError as exc:
                print("End the REST service:")
                print(exc)
            _rm(start_file)
        _rm(properties_file)

    def annotation(self, gene: str, what: str) -> str:
        """ExomiserAnnotation: what = "all" or "score"."""
        if what == "all":
            return self.table.get(gene, "-1.0\t\t\t")
        if what == "score":
            return self.table[gene].split("\t")[0] if gene in self.table else "-1.0"
        print(f'proc ExomiserAnnotation: Bad option value for "what" ({what})')
        return ""
