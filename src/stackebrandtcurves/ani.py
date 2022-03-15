import argparse
import collections
import itertools
import os
import os.path
import random
import re
import shutil
import subprocess
import urllib.error

from .download import get_url
from .parse import parse_fasta, write_fasta



def remove_files(target_dir):
    if os.path.exists(target_dir):
        for filename in os.listdir(target_dir):
            os.remove(os.path.join(target_dir, filename))

class AssemblyPair:
    genome_dir = "temp_genomes"
    ani_dir = "temp_ani"
    cache_dir = "temp_genomes_cache"
    ani_cache = {}

    def __init__(
            self, query, subject, pctid=None,
            query_seqid=None, subject_seqid=None):
        self.query = query
        self.subject = subject
        self._pctid = pctid
        self.ani = None
        self.query_seqid = query_seqid
        self.subject_seqid = subject_seqid

    @property
    def pctid(self):
        return round(float(self._pctid), 1)

    def compute_ani(self):
        ani_key = (self.query.accession, self.subject.accession)
        if ani_key in self.ani_cache:
            print("ANI is cached")
            self.ani = self.ani_cache[ani_key]
            return

        if not os.path.exists(self.cache_dir):
            os.mkdir(self.cache_dir)

        remove_files(self.genome_dir)
        
        query_fname = "{0}.fna.gz".format(self.query.accession)
        query_genome_fp = os.path.join(self.genome_dir, query_fname)
        query_cache_fp = os.path.join(self.cache_dir, query_fname)
        if os.path.exists(query_cache_fp):
            os.rename(query_cache_fp, query_genome_fp)
        else:    
            self.query.download_genome(self.genome_dir, query_fname)
            shutil.copy(query_genome_fp, query_cache_fp)
        subprocess.check_call(["gunzip", "-f", query_genome_fp])

        subject_fname = "{0}.fna.gz".format(self.subject.accession)
        subject_genome_fp = os.path.join(self.genome_dir, subject_fname)
        self.subject.download_genome(self.genome_dir, subject_fname)
        subprocess.check_call(["gunzip", "-f", subject_genome_fp])

        print(
            "Computing ANI for", self.query.accession, "and",
            self.subject.accession)
        ani_fp = "tmp_ani.txt"
        subprocess.check_call([
            "fastANI",
            "-q", self.query_genome_fp,
            "-r", self.reference_genome_fp,
            "-o", self.ani_fp,
        ])

        with open(ani_fp) as f:
            ani_result = next(parse_pairwise_ani(f))
        print("ANI:", ani_result["ani"])
        self.ani_cache[ani_key] = ani_result["ani"]
        self.ani = ani_result["ani"]

    def format_output(self):
        if self.query_seqid and self.subject_seqid:
            return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(
                self.query.accession, self.subject.accession,
                self.query_seqid, self.subject_seqid,
                self.pctid, self.ani,
            )
        return "{0}\t{1}\t{2}\t{3}\n".format(
            self.query.accession, self.subject.accession,
            self.pctid, self.ani,
        )

ANI_FIELDS = [
    "query_fp", "ref_fp", "ani", "fragments_aligned", "fragments_total"]
ANI_TYPES = [str, str, float, int, int]

def parse_pairwise_ani(f):
    for line in f:
        toks = line.strip().split()
        vals = [fcn(tok) for tok, fcn in zip(toks, ANI_TYPES)]
        yield dict(zip(ANI_FIELDS, vals))


def pctid_range(min_pctid):
    assert 50.0 < min_pctid <= 100.0
    current_val = 100.0
    while current_val > min_pctid:
        yield current_val
        current_val = current_val - 0.1



def main_train_soft_threshold(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument(
        "--output-file", type=argparse.FileType("w"),
        default="refseq_pctid_ani.tsv",
        help="Output file",
    )
    p.add_argument(
        "--min_pctid", type=float, default=97.0,
        help="Minimum 16S percent ID",
    )
    p.add_argument(
        "--num-threads", type=int,
        help="Number of threads for 16S percent ID (default: use all CPUs)",
    )
    p.add_argument(
        "--num-ani", type=int, default=100,
        help="Number of genome pairs on which to evaluate ANI",
    )
    p.add_argument(
        "--seed", type=int, default=42,
        help="Random number seed",
    )
    args = p.parse_args()
    args.output_file.write(
        "query_assembly\tsubject_assembly\t"
        "query_seqid\tsubject_seqid\t"
        "pctid\tani\n")
    
    # Set seed for 16S selection
    random.seed(args.seed)
    
    # Load all assemblies
    assemblies = RefseqAssembly.load()

    # 16S database with one random sequence from each assembly
    db = Refseq16SDatabase(
        "refseq_16S_all.fasta",
        "refseq_16S_accessions_all.txt")
    if os.path.exists(db.accession_fp):
        db.load(assemblies)
    else:
        for assembly in assemblies.values():
            db.add_assembly(assembly)
        db.save()

    pctid_vals = list(pctid_range(args.min_pctid)) * args.num_ani

    # Set seed again
    random.seed(args.seed + 1)
    assembly_list = list(assemblies.values())
    for current_pctid in pctid_vals:
        found = False
        while not found:
            # randomly select assembly
            query_assembly = random.choice(assembly_list)
            # randomly select query sequence from assembly
            query_assembly_seqids = db.seqids_by_assembly[query_assembly.accession]
            # next loop if we don't have any sequences for this assembly
            if not query_assembly_seqids:
                continue
            query_seqid = random.choice(query_assembly_seqids)
            assembly_pairs = db.search_one(
                query_seqid, current_pctid, threads=args.num_threads)
            assembly_pairs = list(assembly_pairs)
            if assembly_pairs:
                try:
                    # randomly select one result
                    selected_pair = random.choice(assembly_pairs)
                    selected_pair.compute_ani()
                    args.output_file.write(selected_pair.format_output())
                    args.output_file.flush()
                except Exception as e:
                    print(e)
                else:
                    found = True

