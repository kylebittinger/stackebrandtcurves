import collections
import io
import os
import re
import shutil
import subprocess
import urllib.request


class RefSeq:
    summary_url = (
        "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/"
        "bacteria/assembly_summary.txt"
        )

    def __init__(self, data_dir="refseq_data", max_n=5):
        self.data_dir = data_dir
        self.max_n = max_n
        self.assemblies = {}
        self.seqs = {}
        self.accession_seqids = collections.defaultdict(list)
        self.seqid_accessions = {}

    @property
    def assembly_summary_fp(self):
        return os.path.join(self.data_dir, "assembly_summary.txt")

    def download_summary(self):
        fp = self.assembly_summary_fp
        if not os.path.exists(fp):
            get_url(self.summary_url, fp)
        return fp

    def load(self):
        self.load_assemblies()
        self.load_seqs()

    def load_assemblies(self):
        if not os.path.exists(self.data_dir):
            os.makedirs(self.data_dir)
        self.download_summary()
        with open(self.assembly_summary_fp) as f:
            for assembly in RefseqAssembly.parse(f):
                self.assemblies[assembly.accession] = assembly
        return self.assemblies

    def load_seqs(self):
        if os.path.exists(self.accession_fp):
            self.reload_seqs()
        else:
            self.collect_seqs()
            self.save_seqs()

    def reload_seqs(self):
        with open(self.ssu_fasta_fp) as f:
            for seqid, seq in parse_fasta(f):
                self.seqs[seqid] = seq
        with open(self.accession_fp) as f:
            for seqid, accession in parse_accessions(f):
                self.seqid_accessions[seqid] = accession
                self.accession_seqids[accession].append(seqid)

    def collect_seqs(self):
        for accession in self.assemblies.keys():
            seqs = list(self.get_16S_seqs(accession))
            for seqid, seq in seqs:
                self.accession_seqids[accession].append(seqid)
                self.seqs[seqid] = seq
                self.seqid_accessions[seqid] = accession

    def save_seqs(self):
        with open(self.ssu_fasta_fp, "w") as f:
            for seqid, seq in self.seqs.items():
                f.write(">{0}\n{1}\n".format(seqid, seq))
        with open(self.accession_fp, "w") as f:
            for seqid, accession in self.seqid_accessions.items():
                f.write("{0}\t{1}\n".format(seqid, accession))

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
        get_url(assembly.rna_url, rna_fp + ".gz")
        subprocess.check_call(["gunzip", "-q", rna_fp + ".gz"])
        return rna_fp

    def get_16S_seqs(self, accession):
        rna_fp = self.download_rna(accession)
        with open(rna_fp, "rt") as f:
            for desc, seq in parse_fasta(f):
                if self.is_16S(desc, seq):
                    print(desc)
                    seqid = desc.split()[0]
                    yield seqid, seq

    def is_16S(self, desc, seq):
        if is_full_length_16S(desc):
            if not is_low_quality(seq, self.max_n):
                return True
        return False

    @property
    def ssu_fasta_fp(self):
        return os.path.join(self.data_dir, "refseq_16S.fasta")

    @property
    def accession_fp(self):
        return os.path.join(self.data_dir, "refseq_16S_accessions.txt")


class RefseqAssembly:
    fields = [
        "assembly_accession", "bioproject", "biosample", "wgs_master",
        "refseq_category", "taxid", "species_taxid", "organism_name",
        "infraspecific_name", "isolate", "version_status", "assembly_level",
        "release_type", "genome_rep", "seq_rel_date", "asm_name", "submitter",
        "gbrs_paired_asm", "paired_asm_comp", "ftp_path",
        "excluded_from_refseq", "relation_to_type_material",
    ]

    def __init__(self, assembly_accession, ftp_path, **kwargs):
        self.accession = assembly_accession
        self.ftp_path = ftp_path
        for key, val in kwargs.items():
            setattr(self, key, val)

    @classmethod
    def parse(cls, f):
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#") or (line == ""):
                continue
            toks = line.split("\t")
            vals = dict(zip(cls.fields, toks))
            if vals["ftp_path"] == "na":
                continue
            yield cls(**vals)

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


def parse_fasta(f):
    f = iter(f)
    try:
        desc = next(f).strip()[1:]
    except StopIteration:
        return
    seq = io.StringIO()
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            yield desc, seq.getvalue()
            desc = line[1:]
            seq = io.StringIO()
        else:
            seq.write(line)
    yield desc, seq.getvalue()


def parse_accessions(f):
    for line in f:
        if line.startswith("#"):
            continue
        line = line.strip()
        yield line.split("\t")


def is_low_quality(seq, max_n):
    ns = "N" * max_n
    return ns in seq

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


def get_url(url, fp):
    print("Downloading", url)
    with urllib.request.urlopen(url) as resp, open(fp, 'wb') as f:
        shutil.copyfileobj(resp, f)
    return fp
