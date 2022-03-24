# stackebrandtcurves

Stackebrandt curves show the relationship between average nucleotide
identity (ANI) and 16S rRNA gene similarity for bacterial genomes.

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
