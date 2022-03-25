import os
import re
import subprocess

from .download import get_url
from .parse import parse_fasta, parse_assembly_summary

class RefSeq:
    summary_url = (
        "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/"
        "bacteria/assembly_summary.txt"
        )

    def __init__(self, data_dir="refseq_data"):
        self.data_dir = data_dir
        self.assemblies = {}
        self.seqs = {}
        self.accession_seqids = {}
        self.seqid_accessions = {}

    @property
    def assembly_summary_fp(self):
        return os.path.join(self.data_dir, "assembly_summary.txt")

    def download_summary(self):
        fp = self.assembly_summary_fp
        if not os.path.exists(fp):
            get_url(self.summary_url, fp)
        return fp

    def load_assemblies(self):
        self.download_summary()
        with open(self.assembly_summary_fp) as f:
            for assembly in RefseqAssembly.parse(f):
                self.assemblies[assembly.accession] = assembly
        return self.assemblies

    def load_seqs(self):
        for accession in self.assemblies.keys():
            seqs = list(self.get_16S_seqs(accession))
            self.accession_seqids[accession] = [seqid for seqid, seq in seqs]
            for seqid, seq in seqs:
                self.seqs[seqid] = seq
                self.seqid_accessions[seqid] = accession

    def save_seqs(self):
        with open(self.ssu_fasta_fp, "w") as f:
            for seqid, seq in self.seqs.items():
                f.write(">{0}\n{1}\n".format(seqid, seq))

    def save_filtered_seqs(self, fp, seen):
        with open(fp, "w") as f:
            for seqid, seq in self.seqs.items():
                if seqid not in seen:
                    f.write(">{0}\n{1}\n".format(seqid, seq))

    @property
    def genome_dir(self):
        return os.path.join(self.data_dir, "genome_fasta")

    def genome_fp(self, assembly):
        genome_filename = "{0}_genomic.fna".format(assembly.basename)
        return os.path.join(self.genome_dir, genome_filename)

    def collect_genome(self, accession):
        assembly = self.assemblies[accession]
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

    def download_rna(self, accession):
        assembly = self.assemblies[accession]
        rna_fp = self.rna_fp(assembly)
        if os.path.exists(rna_fp):
            return rna_fp
        if not os.path.exists(self.rna_dir):
            os.makedirs(self.rna_dir)
        print("Downloading 16S seqs for ", assembly.accession)
        get_url(assembly.rna_url, rna_fp + ".gz")
        subprocess.check_call(["gunzip", "-q", rna_fp + ".gz"])
        return rna_fp

    def get_16S_seqs(self, accession):
        rna_fp = self.download_rna(accession)
        with open(rna_fp, "rt") as f:
            for desc, seq in parse_fasta(f):
                if is_full_length_16S(desc):
                    print(desc)
                    seqid = desc.split()[0]
                    yield seqid, seq

    @property
    def ssu_fasta_fp(self):
        return os.path.join(self.data_dir, "refseq_16S.fasta")


class RefseqAssembly:
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


def is_full_length_16S(desc):
    accession, attrs = parse_desc(desc)
    product = attrs.get("product")
    is_16S = product == "16S ribosomal RNA"
    if not is_16S:
        return False
    location = attrs.get("location", "")
    is_full_length = (">" not in location) and ("<" not in location)
    return is_full_length


def parse_desc(desc):
    accession, sep, rest = desc.partition(" ")
    toks = re.findall(r'\[([^\]]+)\]', rest)
    attrs = {}
    for tok in toks:
        attr, sep, val = tok.partition("=")
        attrs[attr] = val
    return accession, attrs
