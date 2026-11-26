"""Port of ``share/tcl/AnnotSV/AnnotSV-config.tcl`` (``configureAnnotSV``).

Options are loaded in this order, later sources winning:

1. built-in defaults,
2. the ``configfile`` (``etc/AnnotSV/configfile``, or a ``configfile`` sitting
   next to the SV input file),
3. ``-option value`` pairs from the command line.

Option names are case-insensitive. Values are stored as strings (or lists of
strings once validated) under the same keys the Tcl ``g_AnnotSV`` array used,
so later ported modules can be checked line by line against the Tcl code.
"""

from __future__ import annotations

import datetime
import glob
import gzip
import os
import re
import shlex
import subprocess
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence

from . import tcl
from .general import is_an_empty_file, lines_from_file


class ConfigExit(Exception):
    """Raised where the Tcl code called ``exit``; carries the exit code."""

    def __init__(self, message: str, code: int = 2) -> None:
        super().__init__(message)
        self.message = message
        self.code = code


BANNER = "#" * 76

DEFAULT_OPTIONS: Dict[str, object] = {
    "annotationsDir": "",
    "annotationMode": "both",
    "bcftools": "bcftools",
    "bedtools": "bedtools",
    "benignAF": "0.01",
    "bracketedAltMaxSize": "10",
    "bracketedAltMode": "ALL",
    "breakpointProximity": "100",
    "candidateGenesFile": "",
    "candidateGenesFiltering": "0",
    "candidateSnvIndelFiles": "",
    "candidateSnvIndelSamples": "",
    "extann": [],  # external gene annotation files, completed later
    "externalGeneFiles": "",
    "genomeBuild": "GRCh38",
    "hpo": "",
    "includeCI": "1",
    "metrics": "us",
    "minTotalNumber": "500",
    "missingGTinSamplesid": "1",
    "outputColHeader": [],  # not settable from the command line
    "outputDir": "",
    "outputFile": "",
    "overlap": "100",
    "overwrite": "1",
    "promoterSize": "500",
    "rankFiltering": "1-5,NA",
    "ranking": "1",  # not settable from the command line
    "reciprocal": "0",
    "REreport": "0",
    "REselect1": "1",
    "REselect2": "1",
    "samplesidBEDcol": "-1",
    "samplesidTSVcol": "-1",  # not settable from the command line
    "snvIndelFiles": "",
    "snvIndelPASS": "0",
    "snvIndelSamples": "",
    "SVinputFile": "",
    "SVinputInfo": "1",
    "SVminSize": "50",
    "svtBEDcol": "-1",
    "svtTSVcol": "-1",  # not settable from the command line
    "tx": "RefSeq",
    "txFile": "",
    "variantconvertDir": "",
    "variantconvertMode": "combined",
    "vcf": "0",
}

DEFAULT_RANKING_CRITERIA: Dict[str, float] = {
    "1A_Loss": 0.00, "1B_Loss": -0.60, "2A_Loss": 1.00, "2B_Loss": 0.00,
    "2C-1_Loss": 0.90, "2C-2_Loss": 0.00, "2D-1_Loss": 0.00, "2D-2_Loss": 0.90,
    "2D-3_Loss": 0.45, "2D-4_Loss": 0.90, "2E-1_Loss": 0.90, "2E-2_Loss": 0.45,
    "2E-3_Loss": 0.30, "2E-4_Loss": 0.20, "2F-4O_Loss": -1.00, "2G_Loss": 0.00,
    "2H_Loss": 0.15, "3A_Loss": 0.00, "3B_Loss": 0.45, "3C_Loss": 0.90,
    "4O_Loss": -1.00, "5F_Loss": 0.00, "5G_Loss": 0.10, "5H_Loss": 0.30,
    "1A_Gain": 0.00, "1B_Gain": -0.60, "2A_Gain": 1.00, "2B_Gain": 0.00,
    "2C_Gain": -1.00, "2D_Gain": -1.00, "2E_Gain": 0.00, "2F_Gain": -1.00,
    "2G_Gain": 0.00, "2H-1_Gain": 0.45, "2H-2_Gain": 0.00, "2I-1_Gain": 0.45,
    "2I-2_Gain": 0.45, "2I-3_Gain": 0.00, "2J_Gain": 0.00, "2K_Gain": 0.45,
    "2L_Gain": 0.00, "3A_Gain": 0.00, "3B_Gain": 0.45, "3C_Gain": 0.90,
    "4O_Gain": 0.00, "5F_Gain": 0.00, "5G_Gain": 0.10, "5H_Gain": 0.15,
}

# Options a user may set (configfile or command line), in Tcl order.
OPTIONS_OK: List[str] = (
    "annotationsDir annotationMode bcftools bedtools benignAF bracketedAltMaxSize "
    "bracketedAltMode breakpointProximity candidateGenesFile candidateGenesFiltering "
    "candidateSnvIndelFiles candidateSnvIndelSamples extann externalGeneFiles genomeBuild "
    "hpo includeCI metrics minTotalNumber missingGTinSamplesid outputDir outputFile overlap "
    "overwrite promoterSize rankFiltering reciprocal REreport REselect1 REselect2 "
    "samplesidBEDcol snvIndelFiles snvIndelPASS snvIndelSamples SVinputFile SVinputInfo "
    "SVminSize svtBEDcol tx txFile variantconvertDir variantconvertMode vcf"
).split()
_OPTIONS_LOWER = {name.lower(): name for name in OPTIONS_OK}

BOOLEAN_OPTIONS = (
    "candidateGenesFiltering includeCI missingGTinSamplesid overwrite reciprocal REreport "
    "REselect1 REselect2 SVinputInfo snvIndelPASS vcf"
).split()

# Columns the user cannot remove from the output (appended if missing).
ESSENTIAL_COLUMNS = (
    "SV_chrom SV_start SV_end SV_length SV_type Annotation_mode Gene_name Gene_Count "
    "RE_gene P_gain_coord P_loss_coord P_snvindel_nb B_gain_coord B_loss_coord "
    "po_P_gain_coord po_P_loss_coord po_B_gain_allG_coord po_B_loss_allG_coord "
    "po_B_loss_someG_coord HI TS GnomAD_pLI LOEUF_bin DDD_HI_percent "
    "PhenoGenius_specificity Exomiser_gene_pheno_score OMIM_morbid Location Location2 "
    "Overlapped_CDS_percent Frameshift Exon_count Samples_ID"
).split()

_GENOME_BUILD_RE = r"^(GRCh37)|(GRCh38)|(CHM13)|(mm39)|(mm9)|(mm10)$"
_HUMAN_RE = r"^(GRCh37)|(GRCh38)|(CHM13)$"
_MOUSE_RE = r"^(mm39)|(mm9)|(mm10)$"


@dataclass
class Config:
    """Replacement for the Tcl ``g_AnnotSV`` and ``g_rankingCriteria`` arrays.

    Index it like the Tcl array: ``config["SVinputFile"]``.
    """

    values: Dict[str, object] = field(default_factory=dict)
    ranking_criteria: Dict[str, object] = field(default_factory=dict)

    def __getitem__(self, key: str):
        return self.values[key]

    def __setitem__(self, key: str, value) -> None:
        self.values[key] = value

    def __contains__(self, key: str) -> bool:
        return key in self.values

    def get(self, key: str, default=None):
        return self.values.get(key, default)


# --------------------------------------------------------------------------
# Install layout (bin/AnnotSV)
# --------------------------------------------------------------------------
def default_install_dir() -> str:
    """Directory holding ``etc/``, ``share/`` and ``commandLineOptions.txt``.

    ``$ANNOTSV`` wins when set; otherwise the repository root this package
    lives in (``<root>/src/annotsv``); otherwise the current directory.
    """
    env = os.environ.get("ANNOTSV")
    if env:
        return os.path.abspath(env.rstrip("/") or "/")
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(os.path.dirname(here))
    if os.path.isdir(os.path.join(root, "etc", "AnnotSV")):
        return root
    return os.getcwd()


def layout(install_dir: str) -> Dict[str, str]:
    return {
        "installDir": install_dir,
        "etcDir": f"{install_dir}/etc/AnnotSV",
        "docDir": f"{install_dir}/share/doc/AnnotSV",
        "bashDir": f"{install_dir}/share/bash/AnnotSV",
    }


# --------------------------------------------------------------------------
# Parsing
# --------------------------------------------------------------------------
def _strip_option_name(name: str) -> str:
    # Tcl: regsub "^-|:[ \t]" (first match only). We also accept "--name"
    # for the Python CLI convenience.
    if name.startswith("--"):
        name = name[1:]
    return re.sub(r"^-|:[ \t]", "", name, count=1)


def _canonical(name: str) -> Optional[str]:
    return _OPTIONS_LOWER.get(name.lower())


def _unknown(name: str, configfile: bool) -> ConfigExit:
    lines = [f'"{name}" option not known.', "For more information on the arguments, please use the -help option"]
    if configfile:
        lines = [BANNER] + lines + [BANNER]
    return ConfigExit("\n".join(lines), 2)


def _pairs(argv: Sequence[str]):
    """Yield ``(name, value)`` pairs like the Tcl ``while {$j < llength}`` loop
    (a trailing option without a value is ignored)."""
    for i in range(0, len(argv) - 1, 2):
        yield _strip_option_name(argv[i]), argv[i + 1]


def find_configfile(config: Config) -> str:
    configfile = f"{config['etcDir']}/configfile"
    input_dir = os.path.dirname(config["SVinputFile"]) or "."
    candidate = f"{input_dir}/configfile"
    if os.path.exists(candidate):
        configfile = candidate
    return configfile


def read_configfile(config: Config, path: str) -> None:
    """Apply options, output columns and ranking weights from a configfile."""
    section = 0  # 0 = options, 1 = output columns, 2 = ranking criteria
    for line in lines_from_file(path):
        if re.search(r"^# AnnotSV Output columns:", line):
            section = 1
        if re.search(r"^# AnnotSV ranking criteria:", line):
            section = 2
        if line.startswith("#") or line == "":
            continue
        if section == 0:
            line = re.sub(r"^-|:[ \t]", "", line)
            items = tcl.split_list(line)
            name = items[0] if items else ""
            value = items[1] if len(items) > 1 else ""
            canonical = _canonical(name)
            if canonical is None:
                raise _unknown(name, configfile=True)
            config[canonical] = value
        elif section == 1:
            config["outputColHeader"].append(re.sub(r"( |\t|\*)+$", "", line))
        else:
            items = tcl.split_list(line)
            if items:
                config.ranking_criteria[items[0]] = items[1] if len(items) > 1 else ""


def apply_arguments(config: Config, argv: Sequence[str]) -> None:
    for name, value in _pairs(argv):
        canonical = _canonical(name)
        if canonical is None:
            raise _unknown(name, configfile=False)
        config[canonical] = value


# --------------------------------------------------------------------------
# Validation
# --------------------------------------------------------------------------
def _bad(*lines: str, code: int = 2) -> ConfigExit:
    return ConfigExit("\n".join([BANNER, *lines, BANNER]), code)


def _resolve_dir(config: Config, key: str, default: str, label: str) -> None:
    if config[key] == "":
        config[key] = default
    else:
        config[key] = re.sub(r"/+$", "", config[key])
        if not os.path.isdir(config[key]):
            raise ConfigExit(
                f"AnnotSV needs in argument an existing path of the {label} directory "
                f'(-{key} = "{config[key]}") - Exit with error.',
                2,
            )
    config[key] = re.sub(r"/+", "/", config[key])


def _vcf_samples(files: List[str], wanted: List[str]):
    """Return (samples found among ``wanted``, all samples) from VCF headers."""
    correct: List[str] = []
    all_samples: List[str] = []
    for path in files:
        opener = gzip.open if re.search(r".gz$", path, re.I) else open
        with opener(path, "rt") as f:
            for line in f:
                line = line.rstrip("\n")
                if not line.startswith("#CHROM"):
                    continue
                ls = line.split("\t")
                all_samples.extend(ls[9:])
                correct.extend(s for s in wanted if s in ls)
    return correct, all_samples


def _resolve_samples(config: Config, files_key: str, samples_key: str) -> None:
    config[samples_key] = tcl.split_chars(config[samples_key], ";|,")
    files = config[files_key]
    if files != "":
        # Tcl: [eval glob -nocomplain $files] (space separated patterns)
        matched: List[str] = []
        for pattern in tcl.split_list(files):
            matched.extend(glob.glob(pattern))
        config[files_key] = matched
    else:
        config[files_key] = []
    if config[files_key]:
        correct, all_samples = _vcf_samples(config[files_key], config[samples_key])
        if not config[samples_key] or not correct:
            config[samples_key] = sorted(set(all_samples))
        else:
            config[samples_key] = sorted(set(correct))
    else:
        config[samples_key] = []


def _check_tools(config: Config, runner: Callable[[List[str]], str]) -> None:
    for tool in ("bedtools", "bcftools"):
        try:
            message = runner(shlex.split(config[tool]) + ["--version"])
        except Exception as exc:  # noqa: BLE001 - any failure means unusable
            raise _bad(f"Bad value for the {tool} option ({config[tool]})", str(exc), "Exit with error.")
        if tool == "bedtools":
            m = re.search(r"([0-9]+\.[0-9]+\.[0-9]+)", message)
            if m:
                major, minor = (int(x) for x in m.group(1).split(".")[:2])
                if major < 2 or (major == 2 and minor < 25):
                    raise _bad(
                        "The minimum bedtools version compatible with AnnotSV is version 2.25.",
                        f"You are using bedtools version {m.group(1)}",
                        "Exit with error.",
                    )


def _run_version(cmd: List[str]) -> str:
    result = subprocess.run(cmd, capture_output=True, text=True)
    output = (result.stdout + result.stderr).strip()
    if result.returncode != 0:
        raise RuntimeError(output or f"exit status {result.returncode}")
    return output


def _is_positive_int_text(value: str) -> bool:
    return re.search(r"[^0-9]", value) is None


def validate(config: Config, check_tools: bool = True,
             tool_runner: Callable[[List[str]], str] = _run_version,
             today: Optional[datetime.date] = None) -> None:
    """Checking of the configuration options (same order as the Tcl code)."""
    c = config
    install_dir = c["installDir"]
    _resolve_dir(c, "annotationsDir", f"{install_dir}/share/AnnotSV", "annotations")
    _resolve_dir(c, "variantconvertDir", f"{install_dir}/share/python3/variantconvert/", "variantconvert")

    # SVinputFile
    if c["SVinputFile"] == "":
        raise ConfigExit(
            "AnnotSV needs in argument the path of your bed or VCF input file (-SVinputFile ...) - Exit with error.", 2)
    if not re.search(r"\.(bed|vcf(.gz)?)$", c["SVinputFile"], re.I):
        raise _bad(f"Bad option value: -SVinputFile = {c['SVinputFile']}",
                   'Extension file should be ".bed" or ".vcf" - Exit with error.')
    if not os.path.exists(c["SVinputFile"]):
        raise _bad(f"Bad value for the SVinputFile option, file does not exist ({c['SVinputFile']}) - Exit with error.")
    if is_an_empty_file(c["SVinputFile"]):
        raise _bad(f"SVinputFile ({c['SVinputFile']} is empty, no SV to annotate - Exit without error.", code=0)

    # bracketedAltMaxSize: positive float
    v = c["bracketedAltMaxSize"]
    if not tcl.string_is_double(v) or tcl.compare(v, "0") < 0:
        raise _bad(f"Bad option value: -bracketedAltMaxSize = {v}", "Should be a positive float, default = 10")

    # benignAF: [0.001-0.1]
    v = c["benignAF"]
    if not tcl.string_is_double(v) or tcl.compare(v, "0.001") < 0 or tcl.compare(v, "0.1") > 0:
        raise _bad(f"Bad option value: -benignAF = {v}", "Should be in the [0.001-0.1] range values, default = 0.01")

    # metrics: us or fr (warning only)
    if c["metrics"] not in ("us", "fr"):
        print(BANNER)
        print(f"Bad option value: -metrics = {c['metrics']}")
        print('Should be "us" or "fr" - Continue with the "us" default.')
        print(BANNER)
        c["metrics"] = "us"

    for key in BOOLEAN_OPTIONS:
        if c[key] not in ("1", "0"):
            raise _bad(f"Bad option value: -{key} = {c[key]}", 'Should be "1" or "0"')

    # hpo: comma, semicolon, pipe or space separated HP:nnnn terms
    correct, display = [], []
    for term in tcl.split_chars(c["hpo"], ",|; "):
        if re.search(r"^HP:[0-9]+$", term):
            correct.append(term)
        else:
            display.append(f"Bad format for the HPO term: {term}. Not used.")
    if display:
        print(BANNER)
        print("\n".join(display))
        print(BANNER)
    c["hpo"] = ",".join(correct)

    # rankFiltering: e.g. "3,4,5", "3-5", "1-5,NA"
    rank = c["rankFiltering"]
    ranks: List[str] = []
    while True:
        m = re.search(r"([1-5])-([1-5])", rank)
        if not m:
            break
        ranks.extend(str(i) for i in range(int(m.group(1)), int(m.group(2)) + 1))
        rank = re.sub(r"[1-5]-[1-5]", "", rank, count=1)
    ranks.extend(i for i in rank.split(",") if i in ("1", "2", "3", "4", "5", "NA"))
    c["rankFiltering"] = sorted(set(ranks))

    # outputDir / outputFile
    if c["outputDir"] == "" and "/" in c["outputFile"]:
        c["outputDir"] = os.path.dirname(c["outputFile"]) or "/"
    if c["outputDir"] != "":
        os.makedirs(c["outputDir"], exist_ok=True)
    else:
        c["outputDir"] = (today or datetime.date.today()).strftime("%Y%m%d") + "_AnnotSV"
        os.makedirs(c["outputDir"], exist_ok=True)
    c["outputDir"] = re.sub(r"/+$", "", c["outputDir"]) or "/"
    if c["outputFile"] != "":
        c["outputFile"] = os.path.basename(c["outputFile"].rstrip("/"))
        if not c["outputFile"].endswith(".tsv"):
            c["outputFile"] += ".tsv"
    else:
        c["outputFile"] = re.sub(r"(\.bed|\.vcf(\.gz)?)$", ".annotated.tsv",
                                 os.path.basename(c["SVinputFile"]), flags=re.I)
    if c["overwrite"] == "0" and os.path.exists(f"{c['outputDir']}/{c['outputFile']}"):
        raise _bad("Bad value for the -outputFile option, file already exists "
                   f"({c['outputDir']}/{c['outputFile']}) - Exit with error.")
    c["outputDir"] = os.path.abspath(c["outputDir"])

    if check_tools:
        _check_tools(c, tool_runner)

    if c["tx"] not in ("RefSeq", "ENSEMBL"):
        raise _bad(f"Bad option value: -tx = {c['tx']}", 'Should be "RefSeq" or "ENSEMBL"')

    if c["txFile"] != "" and not os.path.exists(c["txFile"]):
        raise _bad(f"Bad value for the txFile option, file does not exist ({c['txFile']}) - Exit with error.")

    for key in ("candidateGenesFile", "externalGeneFiles"):
        if c[key] != "" and not os.path.exists(c[key]):
            raise _bad(f"Bad value for the {key} option, file does not exist ({c[key]}) - Exit with error.")

    v = c["overlap"]
    if not _is_positive_int_text(v) or v == "" or int(v) <= 0 or int(v) > 100:
        raise _bad(f"Bad value for the overlap option ({v}), should be an integer ]0,100] - Exit with error.")

    v = c["minTotalNumber"]
    if not _is_positive_int_text(v) or v == "" or int(v) < 100 or int(v) > 1000:
        # The Tcl message names a stale loop variable here; we name the option.
        raise _bad(f"Bad value for the minTotalNumber option ({v}), should be an integer "
                   "comprised between 100 and 1000 - Exit with error.")

    for key in ("breakpointProximity", "promoterSize", "SVminSize"):
        if not _is_positive_int_text(c[key]):
            raise _bad(f"Bad value for the {key} option ({c[key]}), should be a positive integer - Exit with error.")

    if c["annotationMode"] not in ("both", "full", "split"):
        raise _bad(f"Bad option value: -annotationMode = {c['annotationMode']}",
                   'Should be "both", "full" or "split"')

    if c["variantconvertMode"] not in ("combined", "full", "fullsplit"):
        raise _bad(f"Bad option value: -variantconvertMode = {c['variantconvertMode']}",
                   'Should be "combined", "full" or "fullsplit"')

    _resolve_samples(c, "snvIndelFiles", "snvIndelSamples")
    _resolve_samples(c, "candidateSnvIndelFiles", "candidateSnvIndelSamples")
    # TODO(filteredVCF module): the Tcl code then runs
    # `filteredVCFannotation "FULL" ...`, which clears candidateSnvIndelSamples
    # when the candidate VCFs have no GT field. Hook it in when that module is
    # ported.

    # bracketedAltMode: each comma-separated value is used as a regexp against
    # the authorised list, exactly as in Tcl.
    authorised = "ALL BND DEL DUP INS INV"
    for mode in c["bracketedAltMode"].replace(" ", "").split(","):
        try:
            ok = re.search(mode, authorised, re.I) is not None
        except re.error:
            ok = False
        if not ok:
            raise _bad(f"Bad option value: -bracketedAltMode = {c['bracketedAltMode']}",
                       'Should be "ALL", "BND", "DEL", "DUP", "INS" or "INV"')

    if not re.search(_GENOME_BUILD_RE, c["genomeBuild"], re.I):
        raise _bad(f"Bad option value: -genomeBuild = {c['genomeBuild']}",
                   'Should be "GRCh37", "GRCh38", "CHM13", "mm39", "mm9" or "mm10"')
    if re.search(_HUMAN_RE, c["genomeBuild"], re.I):
        c["organism"] = "Human"
    elif re.search(_MOUSE_RE, c["genomeBuild"], re.I):
        c["organism"] = "Mouse"

    c["genesBasedAnn"] = "1"
    for col in ESSENTIAL_COLUMNS:
        if col not in c["outputColHeader"]:
            c["outputColHeader"].append(col)

    # externalGeneFiles is iterated as a Tcl list later on.
    c["externalGeneFiles"] = tcl.split_list(c["externalGeneFiles"])


def configure(argv: Sequence[str], install_dir: Optional[str] = None, check_tools: bool = True,
              tool_runner: Callable[[List[str]], str] = _run_version,
              today: Optional[datetime.date] = None, log=print) -> Config:
    """configureAnnotSV: defaults, then configfile, then ``argv``; then checks.

    Raises :class:`ConfigExit` wherever the Tcl code would ``exit``.
    """
    now = datetime.datetime.now().strftime("%B %d %Y - %H:%M")
    log(f"...downloading the configuration data ({now})")
    log("\t...configuration data by default")
    config = Config(
        values={k: (list(v) if isinstance(v, list) else v) for k, v in DEFAULT_OPTIONS.items()},
        ranking_criteria=dict(DEFAULT_RANKING_CRITERIA),
    )
    config.values.update(layout(install_dir or default_install_dir()))

    for name, value in _pairs(argv):
        if name == "SVinputFile":
            config["SVinputFile"] = value

    configfile = find_configfile(config)
    log(f"\t...configuration data from {configfile}")
    read_configfile(config, configfile)

    log("\t...configuration data given in arguments")
    apply_arguments(config, argv)

    log("\t...checking all these configuration data\n")
    validate(config, check_tools=check_tools, tool_runner=tool_runner, today=today)
    return config


# --------------------------------------------------------------------------
# bin/AnnotSV: external gene-based files and the argument listing
# --------------------------------------------------------------------------
def collect_external_gene_files(config: Config) -> None:
    """Fill ``extann`` / ``geneBasedAnn`` like the tail of ``bin/AnnotSV``."""
    base = f"{config['annotationsDir']}/Annotations_{config['organism']}"
    extann = list(config["extann"])
    gene_based = f"{base}/Gene-based"
    for path in glob.glob(f"{gene_based}/*/*.tsv"):
        if path.endswith("_DGV_samplesInStudies.tsv"):
            continue
        extann.append(path)
    extann.extend(glob.glob(f"{gene_based}/*/*.tsv.gz"))
    user_dir = f"{base}/Users/"
    extann.extend(glob.glob(f"{user_dir}/*.tsv"))
    extann.extend(glob.glob(f"{user_dir}/*.tsv.gz"))
    extann.extend(config["externalGeneFiles"])
    config["extann"] = extann
    config["geneBasedAnn"] = "1" if extann else "0"


_HIDDEN_KEYS = re.compile(
    r"ABCann|Ann|EAann|Exomiser|GHann|MPRAann|PhenoGeniusCli|bashDir|bedFile|docDir|etcDir|"
    r"extann|genes|installDir|outputColHeader|ranking|samplesidTSVcol|svtTSVcol|tclDir|userDir|Version$"
)


def format_value(value) -> str:
    if isinstance(value, list):
        return tcl.format_list(value)
    return str(value)


def arguments_listing(config: Config) -> str:
    """The "AnnotSV has been run with these arguments" block."""
    stars = "\t" + "*" * 42
    lines = ["", "...listing arguments", stars, "\tAnnotSV has been run with these arguments:", stars]
    for key in sorted(config.values):
        if _HIDDEN_KEYS.search(key):
            continue
        value = format_value(config[key])
        if value == "":
            continue
        lines.append(f"\t-{key} {value}")
    lines.append(stars + "\n")
    return "\n".join(lines)
