import argparse
import collections
import itertools
import os
import random

from .assembly import RefseqAssembly
from .ssu import Refseq16SDatabase

def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument(
        "assembly_accession",
        help="Accession of type strain assembly",
    )
    p.add_argument(
        "--output-file",
        help="Output file (default: created from assembly accession)",
    )
    p.add_argument(
        "--min-pctid", type=float, default=90.0,
        help="Minimum 16S percent ID",
    )
    p.add_argument(
        "--max-unique-pctid", type=int, default=100,
        help=(
            "Maximum number of ANI comparisons for each unique value of 16S "
            "percent ID"))
    p.add_argument(
        "--num-threads", type=int,
        help="Number of threads for 16S percent ID (default: use all CPUs)",
    )
    p.add_argument(
        "--seed", type=int, default=42,
        help="Random number seed",
    )
    p.add_argument(
        "--multi-stage-search", action="store_true",
        help="Conduct exhaustive 16S search in several stages",
    )
    p.add_argument(
        "--assembly-summary",
        help="Assembly summary file (default: download from NCBI)",
    )
    args = p.parse_args(argv)

    if args.output_file is None:
        args.output_file = "assembly_{0}_pctid_ani.txt".format(
            args.assembly_accession)
    output_file = open(args.output_file, "w")
    output_file.write(
        "query_assembly\tsubject_assembly\t"
        "query_seqid\tsubject_seqid\t"
        "pctid\tani\tfragments_aligned\tfragments_total\n")

    random.seed(args.seed)

    if args.assembly_summary:
        assembly_summary_fp = args.assembly_summary
    else:
        assembly_summary_fp = "refseq_bacteria_assembly_summary.txt"
        RefseqAssembly.download_summary(assembly_summary_fp)

    with open(assembly_summary_fp, "r") as f:
        assemblies = RefseqAssembly.load(f)

    db = Refseq16SDatabase()
    if os.path.exists(db.accession_fp):
        db.load(assemblies)
    else:
        for assembly in assemblies.values():
            db.add_assembly(assembly)
        db.save()

    query_assembly = assemblies[args.assembly_accession]
    query_assembly_seqids = db.seqids_by_assembly[query_assembly.accession]
    query_seqid = query_assembly_seqids[0]
    query_seq = db.seqs[query_seqid]

    if args.multi_stage_search:
        assembly_pairs = db.exhaustive_search(
            query_seqid, query_seq,
            min_pctid=args.min_pctid,
            threads=args.num_threads)
    else:
        assembly_pairs = db.search_seq(
            query_seqid, query_seq,
            min_pctid=args.min_pctid,
            threads=args.num_threads)
    assembly_pairs = list(assembly_pairs)
    if args.max_unique_pctid:
        pairs_by_pctid = collections.defaultdict(list)
        for pair in assembly_pairs:
            pairs_by_pctid[pair.pctid].append(pair)
        for pctid, pairs in pairs_by_pctid.items():
            if len(pairs) > args.max_unique_pctid:
                pairs_by_pctid[pctid] = random.sample(pairs, k=args.max_unique_pctid)
        assembly_pairs = list(itertools.chain.from_iterable(pairs_by_pctid.values()))
    
    for assembly_pair in assembly_pairs:
        try:
            assembly_pair.compute_ani()
            output_file.write(assembly_pair.format_output())
            output_file.flush()
        except Exception as e:
            print(e)
