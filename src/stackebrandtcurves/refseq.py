import os
import subprocess
import urllib

from .download import get_url
from .parse import parse_fasta

class RefSeq:
    summary_url = (
        "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/"
        "bacteria/assembly_summary.txt"
        )

    def __init__(self, data_dir="refseq_data"):
        self.data_dir = data_dir
        self._16S_seqs = {}
        
    def download_summary(self, fp=None):
        if fp is None:
            fp = os.path.join(self.data_dir, "assembly_summary.txt")
        return get_url(self.summary_url, fp)

    @classmethod
    def parse_summary(cls, f):
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#") or (line == ""):
                continue
            toks = line.split("\t")
            assembly = dict(zip(cls.summary_cols, toks))
            if assembly["ftp_path"] == "na":
                continue
            yield assembly

    @property
    def genome_dir(self):
        return os.path.join(self.data_dir, "genome_fasta")

    def genome_fp(self, assembly):
        genome_filename = "{0}_genomic.fna".format(assembly.basename)
        return os.path.join(self.genome_dir, genome_filename)

    def download_genome(self, assembly):
        genome_fp = self.genome_fp(assembly)
        if os.path.exists(genome_fp):
            return genome_fp
        if not os.path.exists(self.genome_dir):
            os.makedirs(self.genome_dir)
        get_url(assembly.genome_url, genome_fp + ".gz")
        subprocess.check_call(["gunzip", "-q", genome_fp + ".gz"])
        return genome_fp

    @property
    def rna_dir(self):
        return os.path.join(self.data_dir, "rna_fasta")

    def rna_fp(self, assembly):
        rna_filename = "{0}_rna_from_genomic.fna".format(assembly.basename)
        return os.path.join(self.rna_dir, rna_filename)

    def download_rna(self, assembly):
        rna_fp = self.rna_fp(assembly)
        if os.path.exists(rna_fp):
            return rna_fp
        if not os.path.exists(self.rna_dir):
            os.makedirs(self.rna_dir)
        print("Downloading 16S seqs for ", assembly.accession)
        get_url(assembly.rna_url, rna_fp + ".gz")
        subprocess.check_call(["gunzip", "-q", rna_fp + ".gz"])
        return rna_fp

    def get_16S_seqs(self, assembly):
        cached_16S_seqs = self._16S_seqs.get(assembly.accession)
        if cached_16S_seqs is not None:
            return cached_16S_seqs
        rna_fp = self.download_rna(assembly)
        with open(rna_fp, "rt") as f:
            seqs = list(parse_fasta(f))
        res = [(desc, seq) for (desc, seq) in seqs if is_16S(desc)]
        self._16S_seqs[assembly.accession] = res
        return res

def is_16S(desc):
    return "product=16S ribosomal RNA" in desc
