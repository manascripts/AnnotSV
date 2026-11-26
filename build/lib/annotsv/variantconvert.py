"""Port of AnnotSV-variantconvert.tcl (VCF output with variantconvert).

* ``check_variantconvert_configfile`` -> ``checkVariantconvertConfigfile``
* ``run_variantconvert``              -> ``runVariantconvert``

variantconvert itself stays an external tool, run the same way as from Tcl:
``python3 <variantconvertDir>/src/variantconvert convert -i <tsv> -o <vcf> -c <config>``.
"""

from __future__ import annotations

import os
import re

from . import tcl
from .general import content_from_file, lines_from_file, replace_text_in_file, write_text_in_file

_SUFFIX = {"combined": "combined", "full": "full", "fullsplit": "fullsplit"}


def _local_configfile(c, fmt: str) -> str:
    return (f"{c['variantconvertDir']}/src/variantconvert/configs/{c['genomeBuild']}/"
            f"annotsv3_from_{fmt}.{_SUFFIX[c['variantconvertMode']]}.local.json")


def check_variantconvert_configfile(config) -> None:
    """Write the local variantconvert config files (fasta path, mode) and run ``pip install -e .`` once."""
    c = config
    if str(c.get("vcf", "0")) != "1":
        return
    for fmt in ("bed", "vcf"):
        configfile = (f"{c['variantconvertDir']}/src/variantconvert/configs/{c['genomeBuild']}/"
                      f"annotsv3_from_{fmt}.json")
        local = _local_configfile(c, fmt)
        path_dir = c["annotationsDir"] if c["annotationsDir"] != "" else f"{c['installDir']}/share/AnnotSV/"
        fasta = (f"{path_dir}/Annotations_Human/BreakpointsAnnotations/GCcontent/{c['genomeBuild']}/"
                 f"{c['genomeBuild']}_chromFa.fasta")
        dist_path, new_path = r'"path": ".*",', f'"path": "{fasta}",'
        dist_ref, new_ref = r'"##reference=file:.*"', f'"##reference=file:{fasta}"'
        dist_mode = r'"mode": "combined",'
        new_full, new_fullsplit = '"mode": "full",', '"mode": "full\\&split",'

        # A wrong -annotationsDir in an earlier run: the local config file is recomputed
        if os.path.exists(local):
            path_ok = ref_ok = False
            for L in lines_from_file(local):
                m = re.search(r'"path": "(.*)",', L)
                if m and os.path.exists(m.group(1)):
                    path_ok = True
                m = re.search(r'"##reference=file:(.*)"', L)
                if m and os.path.exists(m.group(1)):
                    ref_ok = True
            if not path_ok or not ref_ok:
                os.remove(local)
        if not os.path.exists(local):
            lines = []
            for L in lines_from_file(configfile):
                # Tcl regsub: only the first match, "&" in the new text stands for the match
                L = _regsub(dist_path, L, new_path)
                L = _regsub(dist_ref, L, new_ref)
                if re.search(dist_mode, L):
                    if c["variantconvertMode"] == "full":
                        L = _regsub(dist_mode, L, new_full)
                    elif c["variantconvertMode"] == "fullsplit":
                        L = _regsub(dist_mode, L, new_fullsplit)
                lines.append(L)
            replace_text_in_file("\n".join(lines), local)

    # pip install -e . (also done by the Makefile)
    flag = f"{c['variantconvertDir']}/pipinstall.flag"
    done = False
    try:
        if not os.path.exists(flag) or not re.search("Successfully installed variantconvert",
                                                     content_from_file(flag), re.I):
            for cmd in (["pip3", "install", "-e", "."], ["pip", "install", "-e", "."],
                        ["python", "-m", "pip", "install", "-e", "."]):
                failed, message = tcl.exec_command(cmd, cwd=c["variantconvertDir"])
                if not failed:
                    done = True
                    break
                write_text_in_file(message, flag)
                write_text_in_file("\n##############################################\n\n", flag)
    except OSError:
        pass
    if done:
        write_text_in_file("Done", flag)


def _regsub(pattern: str, line: str, new: str) -> str:
    """Tcl ``regsub`` (first match): "&" is the match, "\\&" a literal "&"."""
    def repl(m):
        return re.sub(r"\\&|&", lambda x: "&" if x.group(0) == "\\&" else m.group(0), new)
    return re.sub(pattern, repl, line, count=1)


def run_variantconvert(config, output_file: str) -> None:
    """runVariantconvert: the VCF output next to the TSV output, and its log file."""
    c = config
    vcf_output = re.sub(r".tsv$", lambda _m: ".vcf", output_file, count=1)
    _, message = tcl.exec_command(["python3", f"{c['variantconvertDir']}/src/variantconvert", "--version"])
    m = re.search(r"variantconvert ([0-9]+\.[0-9]+\.[0-9]+)", message)
    version = f"v{m.group(1)} " if m else ""

    print(f"...creation of the VCF output file: {vcf_output}")
    print(f"   AnnotSV relies on the variantconvert tool {version}(https://github.com/SamuelNicaise/variantconvert)")
    log_file = re.sub(r".vcf$", lambda _m: ".variantconvert.log", vcf_output, count=1)

    fmt = "vcf" if re.search(r"\.vcf(.gz)?$", c["SVinputFile"], re.I) else "bed"
    local = _local_configfile(c, fmt)
    if fmt == "bed" and str(c["svtBEDcol"]) == "-1":
        print('   WARNING: With a "BED" SV input file, the user has to define the -svtBEDcol option.')
        print("            => could not create the VCF output file:")
        print(f"               -svtBEDcol {c['svtBEDcol']}")
        return
    command = f"python3 {c['variantconvertDir']}/src/variantconvert convert -i {output_file} -o {vcf_output} -c {local}"

    _, message = tcl.exec_command(command.split())
    message_reg = message.replace("FutureWarning: Setting an item of incompatible dtype is deprecated and "
                                  "will raise in a future error", "...")
    if re.search("error", message_reg, re.I):
        print("   (a minimal Python 3.8 installation is required, as well as the natsort, panda and pyfaidx "
              "Python modules)")
        print("Error:")
        print(f"   => cf {log_file}")
    replace_text_in_file(f"{command}\n\n{message}", log_file)
