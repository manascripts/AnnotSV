"""Port of the ``bin/AnnotSV`` entry script (the ``annotsv`` console command).

Usage: ``annotsv --sv-input-file in.vcf --option value`` (or ``--option=value``).
The Tcl spellings (``-SVinputFile in.vcf``) still work: option names ignore
case, dashes and underscores, so existing command lines and configfiles run as is.
"""

from __future__ import annotations

import os
import platform
import re
import sys
import time
from typing import List, Optional, Sequence

from . import __version__
from . import regions
from .config import (OPTION_CLI_NAMES, ConfigExit, arguments_listing, collect_external_gene_files,
                     configure, default_install_dir)
from .general import content_from_file
from .closest_genes import ClosestGenes
from .gene_based import GeneBasedAnnotations, check_gene_based_files, check_gene_intolerance_files
from .genes import check_genes_file, genes_annotation
from .regulatory_elements import check_regulatory_elements_files, regulatory_elements_annotation
from .sv_input import read_sv_input
from .variantconvert import check_variantconvert_configfile
from .write import organize_annotation


_HELP_NAMES = dict(OPTION_CLI_NAMES, help="help", version="version")
_HELP_OPTION_RE = re.compile(r"(?<![\w-])-(%s)\b" % "|".join(sorted(_HELP_NAMES, key=len, reverse=True)))


def python_style_help(text: str) -> str:
    """``commandLineOptions.txt`` (shared with the Tcl reference) with the
    Python command and option spellings: ``-SVinputFile:`` -> ``--sv-input-file``."""
    def option_line(m: "re.Match[str]") -> str:
        name = "--" + _HELP_NAMES[m.group(1)]
        width = len(m.group(0))
        return name.ljust(width - 1) + " " if len(name) < width else name + " "

    lines = []
    for line in text.split("\n"):
        line = re.sub(r"\$ANNOTSV/bin/AnnotSV -SVinputFile", "annotsv --sv-input-file", line)
        line = re.sub(r"^-(\w+):[ \t]+", lambda m: option_line(m) if m.group(1) in _HELP_NAMES else m.group(0), line)
        line = _HELP_OPTION_RE.sub(lambda m: "--" + _HELP_NAMES[m.group(1)], line)
        lines.append(line)
    return "\n".join(lines)


def show_help(install_dir: str) -> None:
    """showHelp (AnnotSV-help.tcl): print ``commandLineOptions.txt``."""
    for path in (f"{install_dir}/commandLineOptions.txt",
                 f"{install_dir}/share/doc/AnnotSV/commandLineOptions.txt"):
        if os.path.exists(path):
            print(python_style_help(content_from_file(path)).rstrip("\n"))
            print("\nOption names ignore case, \"-\" and \"_\": the Tcl AnnotSV spellings "
                  "(e.g. -SVinputFile) still work.\n")
            return


def banner(install_dir: str) -> str:
    return "\n".join([
        "",
        "Copyright (C) 2017-current GEOFFROY Veronique",
        "",
        "Please feel free to create a Github issue for any suggestions or bug reports "
        "(https://github.com/lgmgeo/AnnotSV/issues)",
        "",
        f"Python version: {platform.python_version()}",
        "",
        "Application name used:",
        f"{install_dir}\n\n",
    ])


def main(argv: Optional[Sequence[str]] = None, install_dir: Optional[str] = None,
         check_tools: bool = True) -> int:
    args: List[str] = list(sys.argv[1:] if argv is None else argv)
    install_dir = install_dir or default_install_dir()
    joined = " ".join(args)

    print(f"AnnotSV {__version__}")
    if re.search(r"-version|-Version", joined):
        return 0
    print(banner(install_dir))

    if not args:
        print("Arguments are missing see help below\n")
        show_help(install_dir)
        return 0
    # Same (loose) test as the Tcl script: any argument containing "help".
    if re.search("help", joined, re.I) or "-h" in args:
        show_help(install_dir)
        return 0

    try:
        config = configure(args, install_dir=install_dir, check_tools=check_tools)
    except ConfigExit as exc:
        print(exc.message)
        return exc.code

    print(arguments_listing(config))

    try:
        sv = read_sv_input(config)
        check_variantconvert_configfile(config)
        check_genes_file(config)
        check_gene_based_files(config)
        check_gene_intolerance_files(config)
        # Breakpoint annotations, then cytoband
        regions.check_fasta_files(config)
        regions.check_breakpoint_files(config)
        regions.check_cytoband(config)
        check_regulatory_elements_files(config)
        regions.check_tad_files(config)
        regions.check_cosmic_file(config)
        regions.check_users_bed(config)
        # bin/AnnotSV lists the gene-based files once every check* proc has built them
        collect_external_gene_files(config)
        genes = genes_annotation(config)
        regulatory = regulatory_elements_annotation(config, genes.all_genes_overlapped)
        organize_annotation(config, genes, vcf_header=sv.vcf_header, svlen=sv.svlen,
                            closest_genes=ClosestGenes(config), regulatory=regulatory.re,
                            header_file_to_remove=sv.header_file_to_remove,
                            gene_based=GeneBasedAnnotations(config),
                            regions=regions.RegionAnnotations(config),
                            phenogenius=regulatory.phenogenius, exomiser=regulatory.exomiser)
    except ConfigExit as exc:
        print(exc.message)
        return exc.code
    if re.search(r"\.bed$", config["SVinputFile"]):
        header = re.sub(r".bed$", lambda _m: ".header.tsv", config["SVinputFile"], count=1, flags=re.I)
        if os.path.isfile(header):
            os.remove(header)
    print(f"\n...AnnotSV is done with the analysis ({time.strftime('%B %d %Y - %H:%M')})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
