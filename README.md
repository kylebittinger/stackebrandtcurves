# stackebrandtcurves

<!-- Begin badges -->
[![Tests](https://github.com/kylebittinger/stackebrandtcurves/actions/workflows/tests.yml/badge.svg)](https://github.com/kylebittinger/stackebrandtcurves/actions/workflows/tests.yml)
[![CodeCov](https://codecov.io/gh/kylebittinger/stackebrandtcurves/branch/main/graph/badge.svg)](https://codecov.io/gh/kylebittinger/stackebrandtcurves)
[![Super-Linter](https://github.com/kylebittinger/stackebrandtcurves/actions/workflows/linter.yml/badge.svg)](https://github.com/kylebittinger/stackebrandtcurves/actions/workflows/linter.yml)
[![Upload Python Package](https://github.com/kylebittinger/stackebrandtcurves/actions/workflows/python-publish.yml/badge.svg)](https://github.com/kylebittinger/stackebrandtcurves/actions/workflows/python-publish.yml)
<!--  End badges  -->

Stackebrandt curves show the relationship between average nucleotide
identity (ANI) and 16S rRNA gene similarity for bacterial genomes.

## Summary

"Stackebrandt curves" aren't really a thing; I made up the name. The
idea comes from a 1994 paper by Stackebrandt and Goebel[^1] in which
the authors compared 16S rRNA gene similarity and whole-genome
similarity for bacteria. In that paper, the authors used DNA-DNA
hybridization results rather than ANI as a measure of genome
similarity. In a later study, Kim et al.[^2] found that the earlier
results hold true when ANI is used.

To a high degree, 16S gene similarity and full-genome similarity are
correlated. However, Stackebrandt and Goebel observed that 16S
similarity is sometimes high when full-genome similarity is
low. Conversely, full-genome similarity is rarely high if 16S
similarity is low. When plotted with 16S similarity on the y-axis and
full-genome similarity on the x-axis, the relationship looks like a
curve bending along the top of the plot.

Here is an example that follows the classic pattern: *Cutibacterium
acnes* (accession GCF_008728435.1).

![ANI vs. 16S gene similarity for Cutibacterium acnes](/media/cutibacterium_acnes_curve.png)

The purpose of this software is to generate 16S gene similarity and
ANI data for any bacterial genome in NCBI RefSeq. To do this, we
download the list of available genomes, build a 16S gene database,
find other genomes with similar 16S sequences, then compute the ANI
for each one. The program generates a tab-separated output file that
can be used to plot your own Stackebrandt curves.

[^1]: Stackebrandt E and Goebel EM. *International Journal of
Systematic and Environmental Microbiology* **44**, 846 (1994).

[^2]: Kim M, Oh HS, Park SC, Chun J. *International Journal of
Systematic and Environmental Microbiology* **64**, 346 (2014).

## Installation

The Python library and command-line program can be installed using
[pip](https://pypi.org/project/pip/). Besides the python libraries listed 
in the setup file, this program
requires `vsearch` and `fastANI` to be installed. They are both available
through [conda](https://anaconda.org/bioconda/vsearch), which is our
recommended method for installation.

```bash
git clone https://github.com/kylebittinger/stackebrandtcurves.git
cd stackebrandtcurves
conda env create -f stackebrandtcurves_env.yml -n stackebrandtcurves
conda activate stackebrandtcurves
pip install -e .
```

## Usage

The `stackebrandtcurve` program requires one argument, the accesion number
for a bacterial genome assembly in NCBI RefSeq.

```bash
stackebrandtcurve GCF_001688845.2
```

If the program has not been run before, it will automatically download
the files it needs to carry out the computation. The default output file
name is `assembly_GCF_001688845.2_pctid_ani.txt`, but this can be changed
using the `--output-file` option to the program.

Running `stackebrandtcurve --help` will produce a full list of available
options.

## Contributing

We welcome ideas from our users about how to improve this
software. Please open an issue if you have a question or would like to
suggest a feature.
