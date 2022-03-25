from .ssu import limit_results

class StackebrandtApp:
    def __init__(self, db, search_app, ani_app):
        self.db = db
        self.search_app = search_app
        self.ani_app = ani_app
        self.min_pctid = 90.0
        self.max_hits = 100000
        self.threads = None
        self.multi_stage_search = False
        self.max_unique_pctid = 100

    def run(self, query_accession):
        query_seqid = self.db.get_query_seqid(query_accession)
        hits = self.search(query_seqid)

        seqids = [hit["sseqid"] for hit in hits]
        accessions = [self.db.seqid_accessions[s] for s in seqids]
        print(accessions)
        ani_results = self.calculate_ani(query_accession, accessions)

        for hit, ani_result in zip(hits, ani_results):
            seqid = hit["sseqid"]
            accession = self.db.seqid_accessions[seqid]
            yield AppResult(
                query_accession, accession, query_seqid, seqid,
                hit, ani_result)

    def search(self, query_seqid):
        query_seq = self.db.seqs[query_seqid]
        
        if self.multi_stage_search:
            results = self.search_app.exhaustive_search(
                query_seqid, query_seq,
                min_pctid=self.min_pctid, max_hits=self.max_hits,
                threads=self.threads)
        else:
            results = self.search_app.search_seq(
                query_seqid, query_seq,
                min_pctid=self.min_pctid, max_hits=self.max_hits,
                threads=self.threads)

        results = limit_results(results, self.max_unique_pctid)
        hits = [r.hit for r in results]
        return hits

    def calculate_ani(self, query_accession, subject_accessions):
        query_fp = self.db.collect_genome(query_accession)

        unique_subjects = set(subject_accessions)
        subject_fps = {
            self.db.collect_genome(a): a for a in unique_subjects}

        ani_results = self.ani_app.get_ani(query_fp, subject_fps.keys())

        subject_results = {s: None for s in unique_subjects}
        for res in ani_results:
            subject_fp = res["ref_fp"]
            subject = subject_fps[subject_fp]
            subject_results[subject] = res

        return [subject_results[a] for a in subject_accessions]


class AppResult:
    OUTPUT_FIELDS = [
        "query_assembly", "subject_assembly", "query_seqid", "subject_seqid",
        "pctid", "ani", "fragments_aligned", "fragments_total"
    ]
    OUTPUT_TYPES = [str, str, str, str, float, float, int, int]

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
            vals = [fcn(tok) for tok, fcn in zip(toks, cls.OUTPUT_TYPES)]
            yield dict(zip(cls.OUTPUT_FIELDS, vals))

    @classmethod
    def write_header(cls, f):
        f.write("\t".join(cls.OUTPUT_FIELDS))
        f.write("\n")

    def format_output(self):
        pident = self.hit["pident"]
        pident_format = round(float(pident), 2)
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
            self.query.accession, self.subject.accession,
            self.query_seqid, self.subject_seqid,
            pident_format, self.ani_result["ani"],
            self.ani_result["fragments_aligned"],
            self.ani_result["fragments_total"],
        )
