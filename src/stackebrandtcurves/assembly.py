import os
import re
import subprocess
import urllib

from .download import get_url
from .parse import parse_assembly_summary


class RefseqAssembly:
    """A bacterial genome assembly from NCBI RefSeq

    The assembly class holds info on the assebly, and can deliver
    genome or gene sequences for the assembly.
    """
    summary_fp = "refseq_bacteria_assembly_summary.txt"

    summary_url = (
        "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/"
        "bacteria/assembly_summary.txt"
        )
    def __init__(self, assembly_accession, ftp_path, **kwargs):
        self.accession = assembly_accession
        self.ftp_path = ftp_path
        for key, val in kwargs.items():
            setattr(self, key, val)

    @classmethod
    def parse(cls, f):
        for rec in parse_assembly_summary(f):
            yield cls(**rec)

    @property
    def base_url(self):
        return re.sub("^ftp://", "https://", self.ftp_path)

    @property
    def basename(self):
        return os.path.basename(self.ftp_path)

    @property
    def rna_url(self):
        return "{0}/{1}_rna_from_genomic.fna.gz".format(
            self.base_url, self.basename)

    @property
    def genome_url(self):
        return "{0}/{1}_genomic.fna.gz".format(
            self.base_url, self.basename)
