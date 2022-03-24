import collections
import os
import random
import subprocess
import tempfile


class SearchApplication:
    def __init__(self, db, work_dir=None):
        self.db = db
        if work_dir is not None:
            if not os.path.exists(work_dir):
                os.makedirs(work_dir)
            self.work_dir = work_dir
        else:
            self._work_dir_obj = tempfile.TemporaryDirectory()
            self.work_dir = self._work_dir_obj.name

    def get_temp_fp(self, filename):
        return os.path.join(self.work_dir, filename)

    def exhaustive_search(
            self, query_seqid, query_seq, min_pctid=90.0,
            max_hits=10000, threads=None):
        already_found = set()
        already_found.add(query_seqid)

        hits = self.search_seq(
            query_seqid, query_seq, min_pctid=min_pctid,
            max_hits=max_hits, threads=threads)
        for hit in hits:
            already_found.add(hit.subject_seqid)
            yield hit

        for trial in range(10):
            print("Follow-up search, trial {0}".format(trial + 1))
            filtered_subject_fp = self.get_temp_fp("filtered_subject.fasta")
            self.db.save_filtered_seqs(filtered_subject_fp, already_found)
            hits = self.search_seq(
                query_seqid, query_seq, min_pctid=min_pctid,
                max_hits=max_hits, threads=threads,
                subject_fp=filtered_subject_fp)
            n_hits = 0
            for hit in hits:
                n_hits += 1
                already_found.add(hit.subject_seqid)
                yield hit
            if n_hits == 0:
                break
        print("Exhausted 10 search trials")

    def search_seq(
            self, query_seqid, query_seq, min_pctid=90.0,
            max_hits=100000, threads=None, subject_fp=None):
        clear_db = subject_fp is not None
        if subject_fp is None:
            subject_fp = self.db.ssu_fasta_fp

        hits = self.search_once(
            query_seqid, query_seq, subject_fp, min_pctid=min_pctid,
            max_hits=max_hits, threads=threads, clear_db=clear_db)

        for hit in hits:
            query = self.db.seqid_assemblies[hit["qseqid"]]
            subject = self.db.seqid_assemblies[hit["sseqid"]]
            if query.accession != subject.accession:
                yield SearchResult(
                    query, subject, hit["pident"],
                    hit["qseqid"], hit["sseqid"])

    def search_once(
            self, query_seqid, query_seq, subject_fp, min_pctid=90.0,
            max_hits=100000, threads=None, clear_db=False):
        query_fp = self.get_temp_fp("query.fasta")
        with open(query_fp, "w") as f:
            f.write(">{0}\n{1}\n".format(query_seqid, query_seq))
        hits_fp = self.get_temp_fp("hits.txt")

        aligner = PctidAligner(subject_fp)
        aligner.search(
            query_fp, hits_fp, min_pctid=min_pctid, max_hits=max_hits,
            threads=threads)
        if clear_db:
            aligner.clear_db()

        with open(hits_fp) as f:
            for hit in aligner.parse(f):
                hit["vsearch_pident"] = hit["pident"]
                nt_positions = len(hit["qseq"])
                nt_matches = count_matches(hit["qseq"], hit["sseq"])
                hit["pident"] = 100 * (nt_matches / nt_positions)
                yield hit


class PctidAligner:
    field_names = ["qseqid", "sseqid", "pident", "qseq", "sseq"]
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
            max_hits=10000, threads=None):
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
            "query+target+id2+qrow+trow",
            "--maxaccepts", str(max_hits),
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

class SearchResult:
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
        pctid_format = round(float(self.pctid), 2)
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
            self.query.accession, self.subject.accession,
            self.query_seqid, self.subject_seqid,
            pctid_format, self.ani["ani"],
            self.ani["fragments_aligned"], self.ani["fragments_total"],
        )

def limit_results(results, max_results_pctid=None):
    by_pctid = collections.defaultdict(list)
    for result in results:
        by_pctid[result.pctid].append(result)
    for pctid, pctid_results in by_pctid.items():
        if len(pctid_results) > max_results_pctid:
            pctid_results = random.sample(pctid_results, k=max_results_pctid)
        for result in pctid_results:
            yield result

AMBIGUOUS_BASES = {
    "-": "-",
    "T": "T",
    "C": "C",
    "A": "A",
    "G": "G",
    "R": "AG",
    "Y": "TC",
    "M": "CA",
    "K": "TG",
    "S": "CG",
    "W": "TA",
    "H": "TCA",
    "B": "TCG",
    "V": "CAG",
    "D": "TAG",
    "N": "TCAG",
}

def nucleotides_compatible(nt1, nt2):
    if (nt1 == "N") or (nt2 == "N"):
        return True
    a1 = AMBIGUOUS_BASES[nt1]
    a2 = AMBIGUOUS_BASES[nt2]
    return bool(set(a1).intersection(a2))

def nucleotides_match(q, s):
    return (q == s) or nucleotides_compatible(q, s)

def count_matches(qseq, sseq):
    return sum([nucleotides_match(q, s) for q, s in zip(qseq, sseq)])
