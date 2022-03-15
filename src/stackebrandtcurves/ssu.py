import collections
import os
import subprocess

from .ani import AssemblyPair
from .parse import parse_fasta, write_fasta

class Refseq16SDatabase:
    def __init__(
            self, fasta_fp="refseq_16S.fasta",
            accession_fp="refseq_16S_accessions.txt"):
        self.fasta_fp = fasta_fp
        self.accession_fp = accession_fp
        self.seqs = {}
        self.assemblies = {}
        self.seqids_by_assembly = collections.defaultdict(list)

    def add_assembly(self, assembly):
        seqs = list(assembly.ssu_seqs)

        # Avoid writing duplicate sequences for the same genome
        seen = set()
        for desc, seq in seqs:
            if seq not in seen:
                print(desc)
                seqid = desc.split()[0]
                self.assemblies[seqid] = assembly
                self.seqs[seqid] = seq
                self.seqids_by_assembly[assembly.accession].append(seqid)
                seen.add(seq)

    def load(self, assemblies):
        with open(self.accession_fp, "r") as f:
            for line in f:
                toks = line.strip().split()
                seqid = toks[0]
                accession = toks[1]
                assembly = assemblies[accession]
                self.assemblies[seqid] = assembly
                self.seqids_by_assembly[assembly.accession].append(seqid)
        with open(self.fasta_fp, "r") as f:
            for seqid, seq in parse_fasta(f):
                self.seqs[seqid] = seq

    def save(self):
        with open(self.fasta_fp, "w") as f:
            write_fasta(f, self.seqs.items())
        with open(self.accession_fp, "w") as f:
            for seqid, assembly in self.assemblies.items():
                f.write("{0}\t{1}\n".format(seqid, assembly.accession))

    def compute_pctids(self, min_pctid=97.0, threads=None):
        aligner = PctidAligner(self.fasta_fp)
        if not os.path.exists(aligner.hits_fp):
            aligner.search(min_pctid=min_pctid, threads=threads)
        with open(aligner.hits_fp) as f:
            for hit in aligner.parse(f):
                query = self.assemblies[hit["qseqid"]]
                subject = self.assemblies[hit["sseqid"]]
                pctid = hit["pident"]
                yield AssemblyPair(query, subject, pctid)

    def search_one(self, query_seqid, pctid, max_hits=10000, threads=None):
        pctid_str = "{:.1f}".format(pctid)
        print("Searching", query_seqid, "at", pctid_str, "pct identity")
        query_seq = self.seqs[query_seqid]
        query_fp = "temp_query.fasta"
        if os.path.exists(query_fp):
            os.rename(query_fp, "temp_prev_query.fasta")
        query_hits_fp = "temp_query_hits.txt"
        if os.path.exists(query_hits_fp):
            os.rename(query_hits_fp, "temp_prev_query_hits.txt")
        with open(query_fp, "w") as f:
            write_fasta(f, [(query_seqid, query_seq)])
        aligner = PctidAligner(self.fasta_fp)
        aligner.search(
            query_fp, query_hits_fp, min_pctid=pctid,
            threads=threads, max_hits=max_hits)
        with open(query_hits_fp) as f:
            hits = aligner.parse(f)
            for hit in hits:
                if hit["pident"] == pctid_str:
                    query = self.assemblies[hit["qseqid"]]
                    subject = self.assemblies[hit["sseqid"]]
                    pctid = hit["pident"]
                    yield AssemblyPair(
                        query, subject, pctid,
                        hit["qseqid"], hit["sseqid"])

    def search_seq(
            self, query_seqid, query_seq, min_pctid=90.0, max_hits=10000, threads=None):
        query_fp = "temp_query.fasta"
        if os.path.exists(query_fp):
            os.rename(query_fp, "temp_prev_query.fasta")
        query_hits_fp = "temp_query_hits.txt"
        if os.path.exists(query_hits_fp):
            os.rename(query_hits_fp, "temp_prev_query_hits.txt")
        with open(query_fp, "w") as f:
            write_fasta(f, [(query_seqid, query_seq)])
        aligner = PctidAligner(self.fasta_fp)
        aligner.search(
            query_fp, query_hits_fp, min_pctid=min_pctid,
            threads=threads, max_hits=max_hits)
        with open(query_hits_fp) as f:
            hits = aligner.parse(f)
            for hit in hits:
                query = self.assemblies[hit["qseqid"]]
                subject = self.assemblies[hit["sseqid"]]
                pctid = hit["pident"]
                if query.accession != subject.accession:
                    yield AssemblyPair(
                        query, subject, pctid,
                        hit["qseqid"], hit["sseqid"])


class PctidAligner:
    field_names = ["qseqid", "sseqid", "pident"]
    hits_fp = "refseq_16S_hits.txt"

    def __init__(self, fasta_fp):
        self.fasta_fp = fasta_fp

    @property
    def reference_udb_fp(self):
        base_fp, _ = os.path.splitext(self.fasta_fp)
        return base_fp + ".udb"

    def make_reference_udb(self):
        if os.path.exists(self.reference_udb_fp):
            return None
        args = [
            "vsearch",
            "--makeudb_usearch", self.fasta_fp,
            "--output", self.reference_udb_fp,
        ]
        return subprocess.check_call(args)

    def search(
            self, input_fp=None, hits_fp=None, min_pctid=97.0,
            threads=None, max_hits=10000):
        if input_fp is None:
            input_fp = self.fasta_fp
        if hits_fp is None:
            hits_fp = self.hits_fp
        self.make_reference_udb()
        # 97.0 --> 0.97
        min_id = "{:.3f}".format(min_pctid / 100)
        args =[
            "vsearch", "--usearch_global", input_fp,
            "--db", self.reference_udb_fp,
            "--userout", hits_fp,
            "--iddef", "2", "--id", min_id,
            "--userfields",
            "query+target+id2",
        ]
        if max_hits is None:
            args.extend([
                "--maxaccepts", "0", "--maxrejects", "0",
            ])
        else:
            args.extend([
                "--maxaccepts", str(max_hits),
            ])
        if threads is not None:
            args.extend(["--threads", str(threads)])
        print(args)
        subprocess.check_call(args)
        return hits_fp

    def parse(self, f):
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue
            vals = line.split("\t")
            hit = dict(zip(self.field_names, vals))
            if hit["qseqid"] != hit["sseqid"]:
                yield hit
