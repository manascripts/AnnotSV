<p align="center">
    <img src="share/doc/AnnotSV/Images/AnnotSV_logo.png" width="500">

<br />

<div align="center">
    <h1 style="font-weight: bold">An integrated tool for Structural Variations annotation and ranking</h1>
</div>

<br />

# Table of contents

- ## AnnotSV Annotations Engine
    - ### [README](README.AnnotSV_3.5.7.pdf)
    - ### [Getting Started](share/doc/AnnotSV/quickstart.md)
    - ### [Ranking](share/doc/AnnotSV/ranking.md)
    - ### [Annotations](share/doc/AnnotSV/annotations.md)
- ## AnnotSV Webserver
    - ### [Overview](share/doc/AnnotSV/webserver_overview.md)
- ## AnnotSV issues
	<b>Please detail your request when <a href="https://github.com/lgmgeo/AnnotSV/issues" target="_blank">reporting a bug</a>:</b>
	- Give the version of AnnotSV and the command line used
	- Add/join an input file that allow to reproduce the bug
	- Give all the information that could have an impact on the bug reported (modifications made in the installation or in the annotation files)
	- Specify if you use Singularity...
<br />

# Python version

This fork ports AnnotSV from Tcl/shell to Python (package `annotsv`, in `src/annotsv`).
The command line and the output are the same as the Tcl version.

### Install

```bash
git clone https://github.com/manascripts/AnnotSV.git
cd AnnotSV
python3 -m pip install .          # or "pip install -e ." to work on the code
```

This installs two commands:

- `AnnotSV`, which replaces `bin/AnnotSV` (`python3 -m annotsv` works too)
- `AnnotSV-install`, which replaces `bin/INSTALL_annotations.sh` and the data targets of the `Makefile`

`bedtools` (and `bcftools` for `-snvIndelFiles`) must be in your `PATH`, as with the Tcl version. Exomiser also needs `java`.

### Annotations

```bash
AnnotSV-install human                    # make install-human-annotation (includes the Exomiser data)
AnnotSV-install mouse                    # make install-mouse-annotation
AnnotSV-install exomiser-jar             # make install-exomiser-resources
AnnotSV-install variantconvert           # make install-variantconvert (VCF output)
```

By default the annotations go in `<install dir>/share/AnnotSV`, where `AnnotSV` looks for them. The install dir is
`$ANNOTSV` when set, the source tree for an editable install, otherwise the Python environment prefix (where `pip`
puts `etc/AnnotSV/configfile`). For a separate directory (HPC, containers), use the same path for both commands:

```bash
AnnotSV-install --annotations-dir ./AnnotSV_annotations human
AnnotSV -SVinputFile input.vcf -annotationsDir ./AnnotSV_annotations
```

Versions can be changed with `--human-version`, `--exomiser-version`, etc. (`AnnotSV-install --help`).

### Usage

```bash
AnnotSV -SVinputFile input.vcf -outputFile output.annotated.tsv -genomeBuild GRCh38
AnnotSV -help                            # all the options (commandLineOptions.txt)
```

### Tests

```bash
python3 -m pip install -e ".[test]"
python3 -m pytest                        # parity tests run only where tclsh and bedtools are installed
```

<br />

# Collaborative work

Anyone interested in implementing new annotations/features in AnnotSV?

Thanks to the [AnnotSV user community](https://lbgi.fr/AnnotSV/acknowledgments):

    - Bugs could be tackled efficiently

    - New ideas could be investigated faster


We look forward to the opportunity to work together, feel free to fork the page if you want to help :-)


