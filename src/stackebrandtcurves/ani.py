import os
import subprocess
import tempfile

from .parse import parse_pairwise_ani

class AniApplication:
    def __init__(self, db, work_dir=None):
        self.db = db
        self.ani_cache = {}
        if work_dir is not None:
            if not os.path.exists(work_dir):
                os.makedirs(work_dir)
            self.work_dir = work_dir
        else:
            self._work_dir_obj = tempfile.TemporaryDirectory()
            self.work_dir = self._work_dir_obj.name

    def compute_multi_ani(self, search_results):
        query = search_results[0].query
        subjects = list(set(r.subject for r in search_results))

        query_fp = self.db.download_genome(query)
        subject_fps = {self.db.download_genome(s): s for s in subjects}
        reflist_fp = os.path.join(self.work_dir, "ref_list.txt")
        with open(reflist_fp, "w") as f:
            for subject_fp in subject_fps.keys():
                f.write(subject_fp)
                f.write("\n")
        ani_fp = os.path.join(self.work_dir, "ani.txt")

        subprocess.check_call([
            "fastANI",
            "-q", query_fp,
            "--refList", reflist_fp,
            "-o", ani_fp,
        ])

        ani_results = {s: None for s in subjects}
        with open(ani_fp) as f:
            for ani_result in parse_pairwise_ani(f):
                subject = subject_fps[ani_result["ref_fp"]]
                ani_results[subject] = ani_result

        return [ani_results[r.subject] for r in search_results]

    def compute_ani(self, query, subject):
        ani_key = (query.accession, subject.accession)
        ani_result = self.ani_cache.get(ani_key)
        if ani_result is not None:
            print("ANI is cached")
            return ani_result

        query_fp = self.db.download_genome(query)
        subject_fp = self.db.download_genome(subject)

        print(
            "Computing ANI for", query.accession, "and", subject.accession)
        ani_fp = os.path.join(self.work_dir, "temp_ani.txt")
        subprocess.check_call([
            "fastANI",
            "-q", query_fp,
            "-r", subject_fp,
            "-o", ani_fp,
        ])

        with open(ani_fp) as f:
            ani_result = next(parse_pairwise_ani(f))
        print("ANI:", ani_result["ani"])
        self.ani_cache[ani_key] = ani_result
        return ani_result
