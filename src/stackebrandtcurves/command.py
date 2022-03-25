import argparse
import random

from .refseq import RefSeq
from .ssu import SearchApplication, limit_results
from .application import StackebrandtApp, AppResult

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
        "--max-hits", type=int, default=100000,
        help="Maximum number of hits in each search (default: %(default)s)",
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
        "--search-dir",
        help="Directory for search-related files (default: temp directory)",
    )
    p.add_argument(
        "--ani-dir",
        help="Directory for ANI-related files (default: temp directory)",
    )
    p.add_argument(
        "--data-dir", default="refseq_data",
        help="Data directory (default: refseq_data)",
    )
    args = p.parse_args(argv)

    if args.output_file is None:
        args.output_file = "assembly_{0}_pctid_ani.txt".format(
            args.assembly_accession)
    random.seed(args.seed)

    refseq = RefSeq(args.data_dir)
    refseq.load_assemblies()
    refseq.load_seqs()
    refseq.save_seqs()
    search_app = SearchApplication(refseq, work_dir=args.search_dir)

    app = StackebrandtApp(refseq, search_app, args.ani_dir)
    app.min_pctid = args.min_pctid
    app.max_hits = args.max_hits
    app.max_unique_pctid = args.max_unique_pctid
    app.threads = args.num_threads
    app.multi_stage_search = args.multi_stage_search

    results = app.run(args.assembly_accession)

    with open(args.output_file, "w") as f:
        f.write(AppResult.output_header)
        for result in results:
            f.write(result.format_output())
