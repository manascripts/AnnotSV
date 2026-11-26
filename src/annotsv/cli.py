"""Port of the ``bin/AnnotSV`` entry script (argument handling only so far).

The Tcl command line is kept as is: ``annotsv -SVinputFile in.vcf -option value``.
Option names are case-insensitive and ``--option`` is accepted too.
"""

from __future__ import annotations

import os
import platform
import re
import sys
from typing import List, Optional, Sequence

from . import __version__
from .config import ConfigExit, arguments_listing, collect_external_gene_files, configure, default_install_dir
from .general import content_from_file
from .closest_genes import ClosestGenes
from .genes import check_genes_file, genes_annotation
from .sv_input import read_sv_input
from .write import organize_annotation


def show_help(install_dir: str) -> None:
    """showHelp (AnnotSV-help.tcl): print ``commandLineOptions.txt``."""
    for path in (f"{install_dir}/commandLineOptions.txt",
                 f"{install_dir}/share/doc/AnnotSV/commandLineOptions.txt"):
        if os.path.exists(path):
            print(content_from_file(path))
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
    if re.search("help", joined, re.I):
        show_help(install_dir)
        return 0

    try:
        config = configure(args, install_dir=install_dir, check_tools=check_tools)
    except ConfigExit as exc:
        print(exc.message)
        return exc.code

    collect_external_gene_files(config)
    print(arguments_listing(config))

    try:
        sv = read_sv_input(config)
        check_genes_file(config)
        genes = genes_annotation(config)
        organize_annotation(config, genes, vcf_header=sv.vcf_header, svlen=sv.svlen,
                            closest_genes=ClosestGenes(config),
                            header_file_to_remove=sv.header_file_to_remove)
    except ConfigExit as exc:
        print(exc.message)
        return exc.code
    if re.search(r"\.bed$", config["SVinputFile"]):
        header = re.sub(r".bed$", lambda _m: ".header.tsv", config["SVinputFile"], count=1, flags=re.I)
        if os.path.isfile(header):
            os.remove(header)

    # The other annotation sources and the ranking are ported in later modules.
    print("Only the core gene annotation columns are ported to Python so far.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
