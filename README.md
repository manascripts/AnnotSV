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
The output is the same as the Tcl version.

### Install

```bash
git clone https://github.com/manascripts/AnnotSV.git
cd AnnotSV
python3 -m pip install .          # or "pip install -e ." to work on the code
```

This installs three commands:

- `annotsv`, which replaces `bin/AnnotSV` (`python3 -m annotsv` works too)
- `annotsv-install`, which replaces `bin/INSTALL_annotations.sh` and the data targets of the `Makefile`
- `annotsv-liftover`, which replaces `share/tcl/AnnotSV/Scripts/lift_Over_a_BED_file.tcl` (`annotsv-liftover BEDfileToLift outputFile chainFile [liftoverTool]`)

`bedtools` (and `bcftools` for `--snv-indel-files`) must be in your `PATH`, as with the Tcl version. Exomiser also needs `java`, PhenoGenius `poetry` and `git`, and `annotsv-liftover` UCSC's `liftOver`.

The Tcl sources (`bin/AnnotSV`, `share/tcl`) and bash helpers (`share/bash`) stay in the repository as the reference for the parity tests; the Python package does not use them.

### Annotations

```bash
annotsv-install human                    # make install-human-annotation (includes the Exomiser data)
annotsv-install mouse                    # make install-mouse-annotation
annotsv-install exomiser-jar             # make install-exomiser-resources
annotsv-install variantconvert           # make install-variantconvert (VCF output)
```

By default the annotations go in `<install dir>/share/AnnotSV`, where `annotsv` looks for them. The install dir is
`$ANNOTSV` when set, the source tree for an editable install, otherwise the Python environment prefix (where `pip`
puts `etc/AnnotSV/configfile`). For a separate directory (HPC, containers), use the same path for both commands:

```bash
annotsv-install --annotations-dir ./AnnotSV_annotations human
annotsv --sv-input-file input.vcf --annotations-dir ./AnnotSV_annotations
```

Versions can be changed with `--human-version`, `--exomiser-version`, etc. (`annotsv-install --help`).

### Usage

```bash
annotsv --sv-input-file input.vcf --output-file output.annotated.tsv --genome-build GRCh38
annotsv --help                           # all the options
```

Options are written in kebab-case (`--sv-input-file`, `--svt-bed-col`, `--re-report`...), and `--option=value` works too.
Option names ignore case, `-` and `_`, so the Tcl AnnotSV spellings (`-SVinputFile input.vcf -genomeBuild GRCh38`)
and existing configfiles keep working unchanged.

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


