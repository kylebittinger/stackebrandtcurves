from .ani import FastAni
from .search import Vsearch, limit_hits

class StackebrandtApp:
    def __init__(self, db, search_dir=None, ani_dir=None):
        self.db = db
        self.search_app = Vsearch(search_dir)
        self.ani_app = FastAni(ani_dir)
        self.min_pctid = 90.0
        self.max_hits = 100000
        self.threads = None
        self.multi_stage_search = False
        self.max_unique_pctid = 100

    def run(self, query_accession):
        hits = self.search(query_accession)
        if not hits:
            return

        seqids = [hit["sseqid"] for hit in hits]
        accessions = [self.db.seqid_accessions[s] for s in seqids]
        ani_results = self.calculate_ani(query_accession, accessions)

        query_seqid = hits[0]['qseqid']
        for hit, ani_result in zip(hits, ani_results):
            seqid = hit["sseqid"]
            accession = self.db.seqid_accessions[seqid]
            yield AppResult(
                query_accession, accession, query_seqid, seqid,
                hit, ani_result)

    def search(self, query_accession):
        if self.multi_stage_search:
            hits = self.exhaustive_search(query_accession)
        else:
            hits = self.regular_search(query_accession)
        hits = limit_hits(hits, self.max_unique_pctid)
        return list(hits)

    def calculate_ani(self, query_accession, subject_accessions):
        query_fp = self.db.collect_genome(query_accession)

        unique_subjects = set(subject_accessions)
        subject_fps = {
            self.db.collect_genome(a): a for a in unique_subjects}

        ani_results = self.ani_app.run(query_fp, subject_fps.keys())

        subject_results = {s: None for s in unique_subjects}
        for res in ani_results:
            subject_fp = res["ref_fp"]
            subject = subject_fps[subject_fp]
            subject_results[subject] = res

        return [subject_results[a] for a in subject_accessions]

    def regular_search(self, query_accession, subject_fp=None):
        clear_db = subject_fp is not None
        if subject_fp is None:
            subject_fp = self.db.ssu_fasta_fp
        query_seqids = self.db.accession_seqids[query_accession]
        query_seqid = query_seqids[0]
        query_seq = self.db.seqs[query_seqid]
        hits = self.search_app.search_once(
            query_seqid, query_seq, subject_fp, min_pctid=self.min_pctid,
            max_hits=self.max_hits, threads=self.threads, clear_db=clear_db)
        hits = [hit for hit in hits if hit['sseqid'] not in query_seqids]
        return hits

    def exhaustive_search(self, query_accession):
        already_found = set()

        hits = self.regular_search(query_accession)
        for hit in hits:
            already_found.add(hit['sseqid'])
            yield hit

        subject_fp = self.search_app.filtered_fp
        for trial in range(10):
            print("Follow-up search, trial {0}".format(trial + 1))
            self.db.save_filtered_seqs(subject_fp, already_found)
            hits = self.regular_search(query_accession, subject_fp)
            n_hits = 0
            for hit in hits:
                n_hits += 1
                already_found.add(hit['sseqid'])
                yield hit
            if n_hits == 0:
                return
        print("Exhausted 10 search trials")

class AppResult:
    output_fields = [
        "query_assembly", "subject_assembly", "query_seqid", "subject_seqid",
        "pctid", "ani", "fragments_aligned", "fragments_total"
    ]
    output_types = [str, str, str, str, float, float, int, int]
    output_header = "\t".join(output_fields) + "\n"

    def __init__(
            self, query_accession, subject_accession,
            query_seqid, subject_seqid, hit, ani_result):
        self.query_accession = query_accession
        self.subject_accession = subject_accession
        self.query_seqid = query_seqid
        self.subject_seqid = subject_seqid
        self.hit = hit
        self.ani_result = ani_result

    @classmethod
    def parse(cls, f):
        next(f) # skip header
        for line in f:
            toks = line.strip().split("\t")
            vals = [fcn(tok) for tok, fcn in zip(toks, cls.output_types)]
            yield dict(zip(cls.output_fields, vals))

    def format_output(self):
        pident = self.hit["pident"]
        pident_format = round(float(pident), 2)
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
            self.query_accession, self.subject_accession,
            self.query_seqid, self.subject_seqid,
            pident_format, self.ani_result["ani"],
            self.ani_result["fragments_aligned"],
            self.ani_result["fragments_total"],
        )
