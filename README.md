# stackebrandtcurves

Stackebrandt curves show the relationship between average nucleotide
identity (ANI) and 16S rRNA gene similarity for bacterial genomes.

## Summary

"Stackebrandt curves" aren't really a thing; I made up the name. The idea
comes from a 1994 paper by Stackebrandt and Goebel[^1], which may have been
the first to compare 16S rRNA gene similarity with whole-genome
similarity. In that paper, the authors used DNA-DNA hybridization results
rather than ANI between sequenced genomes to compare full-genomes.

To a high degree, 16S gene similarity and full-genome similarity are
correlated. However, Stackebrandt and Goebel observed that 16S
similarity is sometimes high when full-genome similarity is
low. Conversely, full-genome similarity is rarely high if 16S
similarity is low. When plotted with 16S similarity on the y-axis and
full-genome similarity on the x-axis, the relationship looks like a
curve bending along the top of the plot.

The purpose of this software is to generate 16S gene similarity and
ANI data for any bacterial genome in NCBI RefSeq. To do this, we
download the list of available genomes, build a 16S gene database,
find other genomes with similar 16S sequences, then compute the ANI
for each one. The program generates a tab-separated output file that
can be used to plot your own Stackebrandt curves.

[^1]: Stackebrandt E and Goebel EM. *International Journal of
Systematic and Environmental Microbiology* **44**, 846 (1994).

## Installation

The Python library and command-line program can be installed using
[pip](https://pypi.org/project/pip/).

```bash
git clone https://github.com/kylebittinger/stackebrandtcurves.git
cd stackebrandtcurves
pip install -e .
```

Besides the python libraries listed in the setup file, this program
requires `vsearch` and `fastANI` to be installed. They are both available
through [conda](https://anaconda.org/bioconda/vsearch), which is our
recommended method for installation.

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
