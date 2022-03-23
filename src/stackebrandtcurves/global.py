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

def pctid_range(min_pctid):
    assert 50.0 < min_pctid <= 100.0
    current_val = 100.0
    while current_val > min_pctid:
        yield current_val
        current_val = current_val - 0.1
