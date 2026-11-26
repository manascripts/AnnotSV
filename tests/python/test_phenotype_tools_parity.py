"""Compare the Python Exomiser, PhenoGenius and variantconvert integrations with the Tcl procs.

The external tools are mocked (their real output depends on the Java version,
the Exomiser data release, the PhenoGenius model...): fake ``java``,
``poetry``, ``pip3`` and variantconvert scripts are put on the PATH, and both
the Tcl procs and the Python port run against the same mocks.

* Exomiser: ``checkExomiserInstallation`` + ``runExomiser`` against a fake
  REST service (canned JSON answers per NCBI gene ID), then ``g_Exomiser``
  and ``ExomiserAnnotation``.
* PhenoGenius: ``checkPhenoGeniusCli`` + ``runPhenoGeniusCli`` with a fake
  ``poetry`` writing a canned result file, then ``g_PhenoGeniusCli`` and
  ``PhenoGeniusCliAnnotation`` (with the NCBI gene ID fallback).
* variantconvert: ``checkVariantconvertConfigfile`` + ``runVariantconvert``:
  the local config files, the pip flag file, the VCF output and its log.

Runs only where ``tclsh`` is installed. All data is synthetic.
"""
import json
import os
import re
import shutil
import socket
import stat
import subprocess
import sys

import pytest

from test_genes_parity import ROOT
from annotsv import exomiser as exo_mod
from annotsv import helpers as pg_helpers, tcl
from annotsv.config import Config
from annotsv.exomiser import Exomiser, check_exomiser_installation
from annotsv.phenogenius import PhenoGenius, check_phenogenius_cli
from annotsv.ncbi import search_ncbi_gene_id
from annotsv.variantconvert import check_variantconvert_configfile, run_variantconvert

HERE = os.path.join(os.path.dirname(__file__), "tcl_parity")
pytestmark = pytest.mark.skipif(shutil.which("tclsh") is None, reason="tclsh not installed")

FAKE_JAVA = r'''#!/usr/bin/env python3
import json, os, sys, http.server, urllib.parse
port = int([a for a in sys.argv if a.startswith("--server.port=")][0].split("=", 1)[1])
if os.environ.get("FAKE_JAVA_MODE") == "fail":
    print("***************************\nAPPLICATION FAILED TO START\n***************************", flush=True)
    sys.exit(1)
answers = json.load(open(os.environ["FAKE_EXOMISER"]))
class H(http.server.BaseHTTPRequestHandler):
    def do_GET(self):
        q = urllib.parse.parse_qs(urllib.parse.urlparse(self.path).query)
        a = answers.get(q.get("genes", [""])[0], {"status": 200, "body": '{"results": []}'})
        if q.get("phenotypes", [""])[0] != os.environ["FAKE_HPO"]:
            a = {"status": 400, "body": "bad hpo"}
        body = a["body"].encode()
        self.send_response(a["status"])
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)
    def log_message(self, *args):
        pass
print("Started ExomiserApplication", flush=True)
http.server.HTTPServer(("127.0.0.1", port), H).serve_forever()
'''

FAKE_POETRY = r'''#!/usr/bin/env python3
import os, shutil, sys
if sys.argv[1] == "run":
    out = sys.argv[sys.argv.index("--result_file") + 1]
    shutil.copy(os.environ["FAKE_PG_RESULTS"], out)
    print("PhenoGenius done")
'''

FAKE_VARIANTCONVERT = r'''import os, sys
if "--version" in sys.argv:
    print("variantconvert 2.0.1")
    sys.exit(0)
args = sys.argv[1:]
i, o, c = (args[args.index(k) + 1] for k in ("-i", "-o", "-c"))
with open(o, "w") as f:
    f.write("##fileformat=VCFv4.2\n##config=" + os.path.basename(c) + "\n" + open(i).read())
print("converted " + os.path.basename(i))
if os.environ.get("FAKE_VC_ERROR"):
    print("FutureWarning: Setting an item of incompatible dtype is deprecated and will raise in a future error",
          file=sys.stderr)
    if os.environ["FAKE_VC_ERROR"] == "2":
        print("ValueError: fasta not found", file=sys.stderr)
        sys.exit(1)
'''

FAKE_PIP = r'''#!/bin/sh
echo "Obtaining file://$(pwd)"
echo "Successfully installed variantconvert-2.0.1"
'''

NCBI = "\n".join([
    "HGNC ID\tApproved symbol\tPrevious symbols\tAlias symbols\tNCBI Gene ID",
    "HGNC:1\tFGFR2\tCEK3\tBEK\t2263",
    "HGNC:2\tPKD1\t\tPBP\t5310",
    "HGNC:3\tBRCA1\t\t\t672",
    "HGNC:4\tBADID\t\t\t4F-se9[6:7:8]",
    "HGNC:5\tTP53\t\t\t7157",
    "HGNC:6\tNULLS\t\t\t9001",
    "HGNC:7\tHTTERR\t\t\t9002",
    "HGNC:8\tBADJSON\t\t\t9003",
    "HGNC:9\tNOEVID\t\t\t9004",
    "HGNC:10\tEMPTY\t\t\t9005",
    "HGNC:11\tONLYID\t\t\t9006",
])


def _exe(path, text):
    path.write_text(text)
    path.chmod(path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)


def _free_port():
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def _setup(root, port):
    """Annotations, etc and bash dirs + the fake tools dir; returns (env, dirs)."""
    ann = root / "ann"
    ncbi = ann / "Annotations_Human" / "Gene-based" / "NCBIandHGNCgeneID"
    ncbi.mkdir(parents=True)
    (ncbi / "results.txt").write_text(NCBI + "\n")
    bash = root / "bash"
    bash.mkdir()
    _exe(bash / "searchForAFreePortNumber.sh", f"#!/bin/bash\necho {port}\n")
    tools = root / "tools"
    tools.mkdir()
    _exe(tools / "java", FAKE_JAVA.replace("/usr/bin/env python3", sys.executable))
    _exe(tools / "poetry", FAKE_POETRY.replace("/usr/bin/env python3", sys.executable))
    _exe(tools / "pip3", FAKE_PIP)
    etc = root / "etc"
    etc.mkdir()
    out = root / "out"
    out.mkdir()
    install = root / "install"
    (install / "share" / "python3" / "phenogeniuscli" / "PhenoGeniusCli").mkdir(parents=True)
    env = dict(os.environ, PATH=f"{tools}:{os.environ['PATH']}")
    return env, {"annotationsDir": str(ann), "bashDir": str(bash), "etcDir": str(etc), "outputDir": str(out),
                 "installDir": str(install), "organism": "Human"}


def _tcl(what, dump, env, values):
    args = []
    for k, v in values.items():
        args += [k, " ".join(v) if isinstance(v, list) else v]
    r = subprocess.run(["tclsh", os.path.join(HERE, "phenotype_tools.tcl"), ROOT, what, str(dump)] + args,
                       capture_output=True, text=True, env=env)
    assert r.returncode == 0, (r.stdout + r.stderr)[-3000:]
    return dump.read_text()


def _listing(d):
    """File names without their "<date>-<time>_<pid>_" prefix."""
    return sorted(re.sub(r"^[0-9]{8}-[0-9]{6}_[0-9]+_", "", n) for n in os.listdir(d))


# --------------------------------------------------------------------------
# Exomiser
# --------------------------------------------------------------------------

def _evidence(organism, disease, pairs, **extra):
    d = {"model": {"organism": organism, "diseaseTerm": disease}}
    d["bestModelPhenotypeMatches"] = [{"query": {"id": "HP:1", "label": q}, "match": {"id": "HP:2", "label": m}}
                                      for q, m in pairs]
    d.update(extra)
    return d


ANSWERS = {
    "2263": {"status": 200, "body": json.dumps({"results": [{"score": 0.8790527581872061, "phenotypeEvidence": [
        _evidence("HUMAN", "Craniofacial-skeletal-dermatologic dysplasia",
                  [("Brachydactyly", "Finger syndactyly"), ("Broad thumb", "Broad thumb"),
                   ("Craniosynostosis", "Coronal craniosynostosis")]),
        _evidence("MOUSE", "", [("Broad hallux", "abnormal digit"), ("Brachydactyly", "short digits")]),
        _evidence("FISH", "fin defect", [("Zebra", "Alpha")]),
    ]}, {"score": 0.1}]})},
    "5310": {"status": 200, "body": '{"results": [{"score": 1, "phenotypeEvidence": ['
                                    '{"model": {"organism": "HUMAN", "diseaseTerm": 1e2}, '
                                    '"bestModelPhenotypeMatches": [{"query": {"label": "Renal cysts"}, '
                                    '"match": {"id": "HP:3"}}]}, '
                                    '{"model": {"organism": "RAT", "diseaseTerm": "rat disease"}}, '
                                    '{"model": {"diseaseTerm": "no organism"}}, '
                                    '{"model": {"organism": "MOUSE"}, "bestModelPhenotypeMatches": []}]}]}'},
    "672": {"status": 200, "body": json.dumps({"results": [{"score": 0.00004999, "phenotypeEvidence": []}]})},
    "7157": {"status": 200, "body": json.dumps({"results": [{"score": 0.49995, "phenotypeEvidence": [
        _evidence("HUMAN", "Li-Fraumeni syndrome", [("Tumor", "Tumor"), ("Tumor", "Sarcoma")])]}]})},
    "9001": {"status": 200, "body": '{"results": [{"score": 0.5, "phenotypeEvidence": [{"model": '
                                    '{"organism": "HUMAN", "diseaseTerm": null}, "bestModelPhenotypeMatches": '
                                    '[{"query": {"label": null}, "match": null}]}]}]}'},
    "9002": {"status": 500, "body": json.dumps({"results": [{"score": 0.75, "phenotypeEvidence": []}]})},
    "9003": {"status": 200, "body": "<html>not json</html>"},
    "9004": {"status": 200, "body": json.dumps({"results": [{"score": 0.9}]})},
    "9005": {"status": 200, "body": json.dumps({"results": []})},
}

GENES = ["FGFR2", "PKD1", "BRCA1", "BADID", "TP53", "NULLS", "HTTERR", "BADJSON", "NOEVID", "EMPTY", "NOTINHGNC",
         "BEK", "ONLYID"]
HPO = "HP:0001156,HP:0001363,HP:0011304"

EXO_CASES = {  # name: (hpo, java mode, hpo dirs, properties in etc)
    "ok": (HPO, "ok", ["2402", "999", "abc"], False),
    "etcProperties": (HPO, "ok", ["2402"], True),
    "noHpo": ("", "ok", ["2402"], False),
    "noData": (HPO, "ok", [], False),
    "noProperties": (HPO, "ok", ["2402"], None),
    "javaFails": (HPO, "fail", ["2402"], False),
}


def _exomiser_python(values):
    c = Config(dict(values))
    check_exomiser_installation(c)
    out = [f"Exomiser={c['Exomiser']}", f"hpoVersion={c.get('hpoVersion', '')}"]
    e = Exomiser(c)
    if c["Exomiser"] == "1":
        e.run(values["genes"])
    out.append(f"Exomiser={c['Exomiser']}")
    out += [f"{k}\t{e.table[k]}" for k in sorted(e.table)]
    out += [f"{g}\t{e.annotation(g, 'all')}\t{e.annotation(g, 'score')}" for g in values["genes"] + ["NOPE"]]
    return "\n".join(out) + "\n"


@pytest.mark.parametrize("case", sorted(EXO_CASES))
def test_exomiser_matches_tcl(tmp_path, monkeypatch, case):
    hpo, mode, hpo_dirs, etc_props = EXO_CASES[case]
    port = _free_port()
    results = []
    for side in ("t", "p"):
        root = tmp_path / side
        root.mkdir()
        env, values = _setup(root, port)
        exo = os.path.join(values["annotationsDir"], "Annotations_Exomiser")
        os.makedirs(exo)
        for d in hpo_dirs:
            if d == "abc":
                open(os.path.join(exo, d), "w").close()
                continue
            os.makedirs(os.path.join(exo, d))
            if etc_props is False:
                with open(os.path.join(exo, d, "application.properties"), "w") as f:
                    f.write(f"server.port=XXXX\nexomiser.data-directory=YYYY\nexomiser.hpo={d}ZZZZ\nXXXX\n")
        if etc_props:
            with open(os.path.join(values["etcDir"], "application.properties"), "w") as f:
                f.write("port=XXXX & YYYY\\1 ZZZZ\n")
        os.makedirs(os.path.join(values["annotationsDir"], "jar"))
        open(os.path.join(values["annotationsDir"], "jar", "exomiser-rest-prioritiser-14.0.0.jar"), "w").close()
        answers = root / "answers.json"
        answers.write_text(json.dumps(ANSWERS))
        env.update(FAKE_EXOMISER=str(answers), FAKE_HPO=hpo, FAKE_JAVA_MODE=mode)
        values.update(hpo=hpo, genes=GENES)
        dump = root / "dump.txt"
        if side == "t":
            got = _tcl("exomiser", dump, env, values)
        else:
            monkeypatch.setattr(exo_mod, "WAIT_SECONDS", 0.2)
            for k, v in env.items():
                monkeypatch.setenv(k, v)
            # the Tcl side runs the fake searchForAFreePortNumber.sh
            monkeypatch.setattr(exo_mod.helpers, "search_for_a_free_port", lambda: (False, str(port)))
            got = _exomiser_python(values)
            dump.write_text(got)
        results.append((got, _listing(values["outputDir"])))
    (t, t_files), (p, p_files) = results
    assert p.split("\n") == t.split("\n")
    assert p_files == t_files
    if case == "ok":
        assert "FGFR2\t0.8791\tBrachydactyly;Broad thumb;Coronal craniosynostosis;" in t


def test_ncbi_gene_id_quirks(tmp_path):
    _, values = _setup(tmp_path, 1)
    c = Config(values)
    assert search_ncbi_gene_id(c, "BEK") == ""        # "Alias symbol" is not a column of the header
    assert search_ncbi_gene_id(c, "") == "9006"       # ... so "" gets the ID of the last line
    assert search_ncbi_gene_id(c, "BADID") == ""
    assert search_ncbi_gene_id(c, "FGFR2") == "2263"


# --------------------------------------------------------------------------
# PhenoGenius
# --------------------------------------------------------------------------

PG_RESULTS = "\n".join([
    "#gene_id\tgene_symbol\trank\tscore\thpo_implicated\thpo_description_implicated\tphenotype_specificity",
    "5310\tPKD1\t1\t3.2\t[{'HP:0000108': 0.7}, {'HP:0005562': 0.5}]\t[{'Renal corticomedullary cysts': 0.7}, "
    "{'Multiple renal cysts': 0.5}]\tA - the reported phenotype is highly specific",
    "2263\tFGFR2_OLD\t2\t2.9\t[]\t[{\"Broad 'thumb'\": 1.25}, {'Craniosynostosis': 0.25}]\tB - the reported",
    "672\tBRCA1\t3\t0.4\t[]\t[]\tD - not specific",
    "7157\tNOTLISTED\t4\t0.3\t[]\t[{'Sarcoma': 0.3}]\tC - x",
    "9006\tSOMETHING\t5\t0.2\t[]\t[{'Only by id': 0.2}]\t",
    "",
])

PG_CASES = {"ok": (HPO, "0"), "noHpo": ("", "0"), "checkFails": (HPO, "1")}


def _phenogenius_python(values):
    c = Config(dict(values))
    check_phenogenius_cli(c)
    out = [f"PhenoGeniusCli={c['PhenoGeniusCli']}"]
    p = PhenoGenius(c)
    if c["PhenoGeniusCli"] == "1":
        p.run(values["genes"])
    out.append(f"PhenoGeniusCli={c['PhenoGeniusCli']}")
    out += [f"{k}\t{p.table[k]}" for k in sorted(p.table)]
    out += [f"{g}\t{p.annotation(g, 'all')}\t{p.annotation(g, 'specificity')}" for g in values["genes"] + ["NOPE"]]
    return "\n".join(out) + "\n"


@pytest.mark.parametrize("case", sorted(PG_CASES))
def test_phenogenius_matches_tcl(tmp_path, monkeypatch, case):
    hpo, check_exit = PG_CASES[case]
    results = []
    for side in ("t", "p"):
        root = tmp_path / side
        root.mkdir()
        env, values = _setup(root, 1)
        _exe(root / "bash" / "checkPhenoGeniusCliInstall.sh", f"#!/bin/bash\necho checked $1\nexit {check_exit}\n")
        res = root / "pg.tsv"
        res.write_text(PG_RESULTS)
        env.update(FAKE_PG_RESULTS=str(res))
        values.update(hpo=hpo, genes=["PKD1", "FGFR2", "BRCA1", "TP53", "ONLYID", "NOTINHGNC"])
        dump = root / "dump.txt"
        if side == "t":
            got = _tcl("phenogenius", dump, env, values)
        else:
            for k, v in env.items():
                monkeypatch.setenv(k, v)
            # the Tcl side runs the fake checkPhenoGeniusCliInstall.sh
            fake = os.path.join(values["bashDir"], "checkPhenoGeniusCliInstall.sh")
            monkeypatch.setattr(pg_helpers, "check_phenogenius_cli_install",
                                lambda share, fake=fake: tcl.exec_command([fake, share]))
            got = _phenogenius_python(values)
            dump.write_text(got)
        logs = os.listdir(os.path.join(values["installDir"], "share", "python3", "phenogeniuscli"))
        results.append((got, _listing(values["outputDir"]), len(logs)))
    assert results[1][0].split("\n") == results[0][0].split("\n")
    assert results[1][1:] == results[0][1:]
    if case == "ok":
        assert "PKD1\t3.2\tRenal corticomedullary cysts, Multiple renal cysts\tA" in results[0][0]


# --------------------------------------------------------------------------
# variantconvert
# --------------------------------------------------------------------------

def _vc_dir(root, flag):
    vc = root / "variantconvert"
    (vc / "src" / "variantconvert").mkdir(parents=True)
    (vc / "src" / "variantconvert" / "__main__.py").write_text(FAKE_VARIANTCONVERT)
    for build in ("GRCh37", "GRCh38"):
        d = vc / "src" / "variantconvert" / "configs" / build
        d.mkdir(parents=True)
        for fmt in ("bed", "vcf"):
            (d / f"annotsv3_from_{fmt}.json").write_text("\n".join([
                "{", '    "GENOME": {', '        "assembly": "",', '        "path": "/distributed/path.fa",',
                '        "vcf_header": ["##contig", "##reference=file:/distributed/path.fa"]', "    },",
                f'    "GENERAL": {{"origin": "{fmt}",', '        "mode": "combined",', "    }", "}"]) + "\n")
    if flag is not None:
        (vc / "pipinstall.flag").write_text(flag)
    return vc


VC_CASES = {  # name: (input, svtBEDcol, mode, build, pip flag, stale local config, variantconvert error)
    "bed-combined": ("in.bed", "4", "combined", "GRCh38", "Successfully installed variantconvert-2.0\n", False, ""),
    "vcf-full": ("in.vcf.gz", "-1", "full", "GRCh38", "Successfully installed variantconvert\n", False, "1"),
    "vcf-fullsplit": ("in.VCF", "-1", "fullsplit", "GRCh37", "Successfully installed variantconvert\n", True, "2"),
    "bed-noSvt": ("in.bed", "-1", "combined", "GRCh38", "Successfully installed variantconvert\n", False, ""),
    "pipInstall": ("in.bed", "5", "full", "GRCh38", None, False, ""),
}


@pytest.mark.parametrize("case", sorted(VC_CASES))
def test_variantconvert_matches_tcl(tmp_path, monkeypatch, case):
    inp, svt, mode, build, flag, stale, error = VC_CASES[case]
    results = []
    for side in ("t", "p"):
        root = tmp_path / side
        root.mkdir()
        env, values = _setup(root, 1)
        vc = _vc_dir(root, flag)
        fasta = os.path.join(values["annotationsDir"], "Annotations_Human", "BreakpointsAnnotations", "GCcontent",
                             build)
        os.makedirs(fasta)
        open(os.path.join(fasta, f"{build}_chromFa.fasta"), "w").close()
        if stale:
            for fmt in ("bed", "vcf"):
                (vc / "src" / "variantconvert" / "configs" / build / f"annotsv3_from_{fmt}.fullsplit.local.json") \
                    .write_text('"path": "/old/path.fa",\n"##reference=file:/old/path.fa"\n')
        tsv = os.path.join(values["outputDir"], "x.annotated.tsv")
        with open(tsv, "w") as f:
            f.write("AnnotSV_ID\tSV_chrom\n1_1_100_DEL_1\t1\n")
        env.update(FAKE_VC_ERROR=error)
        values.update(vcf="1", variantconvertDir=str(vc), variantconvertMode=mode, genomeBuild=build,
                      SVinputFile=f"/data/{inp}", svtBEDcol=svt, outputFile="x.annotated.tsv")
        if side == "t":
            r = subprocess.run(["tclsh", os.path.join(HERE, "phenotype_tools.tcl"), ROOT, "variantconvert",
                                str(root / "dump.txt")] + [x for kv in values.items() for x in kv],
                               capture_output=True, text=True, env=env, cwd=str(root))
            assert r.returncode == 0, (r.stdout + r.stderr)[-3000:]
            stdout = r.stdout
        else:
            for k, v in env.items():
                monkeypatch.setenv(k, v)
            monkeypatch.chdir(root)
            c = Config(dict(values))
            import io
            import contextlib
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                check_variantconvert_configfile(c)
                run_variantconvert(c, tsv)
            stdout = buf.getvalue()
        files = {}
        for d in (values["outputDir"], str(vc)):
            for dirpath, _, names in os.walk(d):
                for n in names:
                    path = os.path.join(dirpath, n)
                    if n == "__main__.py" or "__pycache__" in path:
                        continue
                    with open(path) as f:
                        text = f.read()
                    files[os.path.relpath(path, str(root))] = text.replace(str(root), "ROOT")
        results.append((files, stdout.replace(str(root), "ROOT")))
    (t, t_out), (p, p_out) = results
    assert sorted(p) == sorted(t)
    for name in t:
        assert p[name].split("\n") == t[name].split("\n"), name
    assert p_out.split("\n") == t_out.split("\n")
    if case == "bed-combined":
        assert "out/x.annotated.vcf" in t
