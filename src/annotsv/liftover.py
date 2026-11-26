"""Port of ``share/tcl/AnnotSV/Scripts/lift_Over_a_BED_file.tcl`` (``AnnotSV-liftover``).

When performing a liftover on a BED file with more than 3 columns (with
annotation columns), it may fail if some annotation values are empty. This:

* performs a liftover on the first three columns (chrom, start, end)
* keeps the associated annotation columns (even if empty)
* adds the source of the chain file (e.g. "b38-GCA_009914755.4.over.chain" =>
  source = "b38") in front of annotations that are genomic positions
  (e.g. 3:12345-45678 => b38_3:12345-45678)

Usage: ``AnnotSV-liftover BEDfileToLift outputFile chainFile [liftoverTool]``
(``liftoverTool`` defaults to ``liftOver``). The temporary files are written in
the current directory, as with the Tcl script.
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
from typing import Dict, List, Optional

from . import tcl
from .general import lines_from_file, write_text_in_file

_POSITION = re.compile(r"^[0-9XYchrMT]+:[0-9]+-[0-9]+")


def chain_source(chain_file: str) -> Optional[str]:
    """``regexp "^(.*?)(-|To|_)" [file tail $chainFile] match source`` without quotes."""
    m = re.match(r"(.*?)(-|To|_)", os.path.basename(chain_file))
    return re.sub("\"|'", "", m.group(1)) if m else None


def _rm(path: str) -> None:
    if os.path.lexists(path):
        os.remove(path)


def lift_over_bed(bed_file: str, output_file: str, chain_file: str, liftover_tool: str = "liftOver") -> int:
    pid = os.getpid()
    to_lift, lifted, unmapped = f"toLift.{pid}.tmp.bed", f"lifted.{pid}.tmp.bed", f"unmapped.{pid}.bed"

    source = chain_source(chain_file)
    if source is None:
        print("...WARNING: No source extracted!")
    else:
        print(f"...genome source: {source}")

    # Parsing of the BED file
    n = 0
    to_write: List[str] = []
    annotations: Dict[str, List[str]] = {}
    print(f"...reading {os.path.basename(bed_file)}")
    print(f"...writing {to_lift}")
    for line in lines_from_file(bed_file):
        if line.startswith("#") or line == "":
            continue
        ls = line.split("\t")
        n += 1
        # BED file with only the chrom, start and end values
        to_write.append("\t".join(ls[0:3]) + f"\t{n}")
        if n % 80000 == 0:
            write_text_in_file("\n".join(to_write), to_lift)
            to_write = []
        # Memorize the annotation columns
        cols = []
        for v in ls[3:]:
            if _POSITION.match(v):
                if source is None:  # the Tcl script stops on its undefined "source" variable
                    print("can't read \"source\": no such variable", file=sys.stderr)
                    return 1
                cols.append(f"{source}_{v}")
            else:
                cols.append(v)
        annotations[str(n)] = cols
    write_text_in_file("\n".join(to_write), to_lift)

    # Liftover
    command = [liftover_tool, to_lift, chain_file, lifted, unmapped]
    print("...running:")
    print(" ".join(command))
    failed, message = tcl.exec_command(command)
    if failed:
        print(message)
        if re.search("error", message, re.IGNORECASE):
            for f in (to_lift, lifted, unmapped):
                _rm(f)
            return 0

    # Check and display different counts
    def count(path, pred):
        try:
            with open(path, "rb") as f:
                return sum(1 for line in f if pred(line))
        except OSError:
            return 0
    n_lifted = count(lifted, lambda line: line.endswith(b"\n"))  # wc -l
    n_unmapped = count(unmapped, lambda line: line.startswith(b"#"))  # grep -c "^#"
    print(f"...{n} genomic positions to lift")
    print(f"...{n_lifted} genomic positions lifted")
    print(f"...{n_unmapped} positions unmapped")
    if n - n_lifted != n_unmapped:
        print("\tplease check these counts")

    # Writing the lifted output file
    to_write = []
    tmp = f"{output_file}.tmp"
    for line in lines_from_file(lifted):
        if line == "":  # Tcl: an empty file has no lines
            continue
        ls = line.split("\t")
        new = ls[0:3]
        id_ = ls[3] if len(ls) > 3 else ""
        new += annotations[id_]
        to_write.append("\t".join(new))
        if id_.isdigit() and int(id_) % 80000 == 0:
            write_text_in_file("\n".join(to_write), tmp)
            to_write = []
    write_text_in_file("\n".join(to_write), tmp)

    # Sort the lifted output file
    with open(output_file, "w") as out:
        subprocess.run(["sort", "-k1,1", "-k2,2n", tmp], stdout=out, check=True)

    # Clean
    print("...cleaning")
    for f in (to_lift, lifted, unmapped, tmp):
        _rm(f)
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) > 1 and os.path.lexists(argv[1]):
        os.remove(argv[1])  # done before the usage check, as in the Tcl script
    if len(argv) < 3:
        print("USAGE:\nAnnotSV-liftover BEDfileToLift outputFile chainFile {liftoverTool}")
        return 0
    tool = argv[3] if len(argv) > 3 and argv[3] != "" else "liftOver"
    return lift_over_bed(argv[0], argv[1], argv[2], tool)


if __name__ == "__main__":
    sys.exit(main())
