import collections
import os
import subprocess

from .parse import parse_fasta, write_fasta

class Refseq16SDatabase:
    def __init__(self, db):
        self.db = db
        self.data_dir = db.data_dir
        self.seqs = {}
        self.assemblies = {}
        self.seqids_by_assembly = collections.defaultdict(list)

    @property
    def fasta_fp(self):
        return os.path.join(self.data_dir, "refseq_16S.fasta")

    @property
    def accession_fp(self):
        return os.path.join(self.data_dir, "refseq_16S_accessions.txt")

    def add_assembly(self, assembly):
        seqs = self.db.get_16S_seqs(assembly)

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

    def load(self, assembly_dict):
        with open(self.accession_fp, "r") as f:
            for line in f:
                toks = line.strip().split()
                seqid = toks[0]
                accession = toks[1]
                assembly = assembly_dict[accession]
                self.assemblies[seqid] = assembly
                self.seqids_by_assembly[assembly.accession].append(seqid)
        with open(self.fasta_fp, "r") as f:
            for seqid, seq in parse_fasta(f):
                self.seqs[seqid] = seq

    def save(self):
        if not os.path.exists(self.data_dir):
            os.makedirs(self.data_dir)
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

    def search_one(self, query_seqid, pctid, threads=None):
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
            threads=threads)
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

    def get_temp_fp(self, filename):
        search_dir = os.path.join(self.data_dir, "temp_search")
        if not os.path.exists(search_dir):
            os.makedirs(search_dir)
        return os.path.join(search_dir, filename)

    def exhaustive_search(
            self, query_seqid, query_seq, min_pctid=90.0, threads=None):
        already_found = set()
        already_found.add(query_seqid)

        hits = self.search_seq(query_seqid, query_seq, min_pctid, threads)
        for hit in hits:
            already_found.add(hit.subject_seqid)
            yield hit

        for trial in range(10):
            print("Follow-up search, trial {0}".format(trial + 1))
            subject_fp = self.get_temp_fp("subject.fasta")
            seqs_to_write = (
                (seq_id, seq) for seq_id, seq in self.seqs.items()
                if seq_id not in already_found)
            with open(subject_fp, "w") as f:
                write_fasta(f, seqs_to_write)
            hits = self.search_seq(
                query_seqid, query_seq, min_pctid, threads, subject_fp)
            n_hits = 0
            for hit in hits:
                n_hits += 1
                already_found.add(hit.subject_seqid)
                yield hit
            if n_hits == 0:
                break
        print("Exhausted 10 search trials")

    # Used by the main command
    def search_seq(
            self, query_seqid, query_seq, min_pctid=90.0, threads=None,
            subject_fp=None):

        if subject_fp is None:
            subject_fp = self.fasta_fp

        query_fp = self.get_temp_fp("query.fasta")
        previous_query_fp = self.get_temp_fp("previous_query.fasta")
        if os.path.exists(query_fp):
            os.rename(query_fp, previous_query_fp)

        hits_fp = self.get_temp_fp("hits.txt")
        previous_hits_fp = self.get_temp_fp("previous_hits.txt")
        if os.path.exists(hits_fp):
            os.rename(hits_fp, previous_hits_fp)

        with open(query_fp, "w") as f:
            write_fasta(f, [(query_seqid, query_seq)])
        aligner = PctidAligner(subject_fp)
        aligner.search(
            query_fp, hits_fp, min_pctid=min_pctid,
            threads=threads)

        with open(hits_fp) as f:
            hits = aligner.parse(f)
            for hit in hits:
                query = self.assemblies[hit["qseqid"]]
                subject = self.assemblies[hit["sseqid"]]
                pctid = hit["pident"]
                if query.accession != subject.accession:
                    yield AssemblyPair(
                        query, subject, pctid,
                        hit["qseqid"], hit["sseqid"])
        if subject_fp is not None:
            aligner.clear_db()


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

    def clear_db(self):
        os.remove(self.reference_udb_fp)

    def search(
            self, input_fp=None, hits_fp=None, min_pctid=97.0,
            threads=None):
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
            "--maxaccepts", str(10000),
        ]
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

class AssemblyPair:
    def __init__(
            self, query, subject, pctid=None,
            query_seqid=None, subject_seqid=None):
        self.query = query
        self.subject = subject
        self.pctid = pctid
        self.ani = None
        self.query_seqid = query_seqid
        self.subject_seqid = subject_seqid

    def format_output(self):
        pctid_format = round(float(self.pctid), 1)
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
            self.query.accession, self.subject.accession,
            self.query_seqid, self.subject_seqid,
            pctid_format, self.ani["ani"],
            self.ani["fragments_aligned"], self.ani["fragments_total"],
        )
