import os
import re
import subprocess
import urllib

from .download import get_url
from .parse import parse_fasta

class RefseqAssembly:
    """A bacterial genome assembly from NCBI RefSeq

    The assembly class holds info on the assebly, and can deliver
    genome or gene sequences for the assembly.
    """
    summary_url = (
        "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/"
        "bacteria/assembly_summary.txt"
        )
    summary_fp = "refseq_bacteria_assembly_summary.txt"
    genome_dir = "genome_fasta"
    rna_dir = "rna_fasta"
    summary_cols = [
        "assembly_accession", "bioproject", "biosample", "wgs_master",
        "refseq_category", "taxid", "species_taxid", "organism_name",
        "infraspecific_name", "isolate", "version_status",
        "assembly_level", "release_type", "genome_rep", "seq_rel_date",
        "asm_name", "submitter", "gbrs_paired_asm", "paired_asm_comp",
        "ftp_path", "excluded_from_refseq", "relation_to_type_material"
    ]

    def __init__(self, assembly_accession, ftp_path, **kwargs):
        self.accession = assembly_accession
        self.ftp_path = ftp_path
        for key, val in kwargs.items():
            setattr(self, key, val)
        self._ssu_seqs = None
            
    @classmethod
    def parse_summary(cls, f):
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#") or (line == ""):
                continue
            toks = line.split("\t")
            vals = dict(zip(cls.summary_cols, toks))
            if vals["ftp_path"] == "na":
                continue
            yield cls(**vals)

    @property
    def ssu_seqs(self):
        if self._ssu_seqs is not None:
            return self._ssu_seqs
        if not os.path.exists(self.rna_fp):
            try:
                self.download_rna()
            except urllib.error.HTTPError as e:
                print(self.accession)
                print(e)
                return []
        with open(self.rna_fp, "rt") as f:
            seqs = list(parse_fasta(f))
        res = [(desc, seq) for (desc, seq) in seqs if is_16S(desc)]
        self._ssu_seqs = res
        return res

    @classmethod
    def load(cls, summary_fp=None):
        if summary_fp is None:
            summary_fp = cls.summary_fp
        if not os.path.exists(summary_fp):
            get_url(cls.summary_url, summary_fp)
        with open(summary_fp, "r") as f:
            return {a.accession: a for a in cls.parse_summary(f)}

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
    def rna_fp(self):
        rna_filename = "{0}_rna_from_genomic.fna".format(self.basename)
        return os.path.join(self.rna_dir, rna_filename)

    def download_rna(self):
        if os.path.exists(self.rna_fp):
            return
        if not os.path.exists(self.rna_dir):
            os.mkdir(self.rna_dir)
        print("Downloading 16S seqs for ", self.accession)
        get_url(self.rna_url, self.rna_fp + ".gz")
        subprocess.check_call(["gunzip", "-q", self.rna_fp + ".gz"])

    @property
    def genome_url(self):
        return "{0}/{1}_genomic.fna.gz".format(
            self.base_url, self.basename)

    @property
    def genome_fp(self):
        genome_filename = "{0}_genomic.fna.gz".format(
            self.basename)
        return os.path.join(self.genome_dir, genome_filename)

    def download_genome(self, genome_dir=None, genome_fname=None):
        if genome_dir is None:
            genome_dir = self.genome_dir
        if genome_fname is None:
            genome_fp = self.genome_fp
        else:
            genome_fp = os.path.join(genome_dir, genome_fname)
        if os.path.exists(genome_fp):
            return
        if not os.path.exists(genome_dir):
            os.mkdir(genome_dir)
        get_url(self.genome_url, genome_fp)


def is_16S(desc):
    return "product=16S ribosomal RNA" in desc

