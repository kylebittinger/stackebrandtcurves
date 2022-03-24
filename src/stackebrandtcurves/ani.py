import os
import subprocess
import tempfile

from .parse import parse_ani

class AniApplication:
    def __init__(self, db, work_dir=None):
        self.db = db
        if work_dir is not None:
            if not os.path.exists(work_dir):
                os.makedirs(work_dir)
            self.work_dir = work_dir
        else:
            self._work_dir_obj = tempfile.TemporaryDirectory()
            self.work_dir = self._work_dir_obj.name

    def get_ani(self, query_genome_fp, subject_genome_fps):
        reflist_fp = os.path.join(self.work_dir, "ref_list.txt")
        with open(reflist_fp, "w") as f:
            for subject_fp in subject_genome_fps:
                f.write(subject_fp)
                f.write("\n")
        ani_fp = os.path.join(self.work_dir, "ani.txt")

        subprocess.check_call([
            "fastANI",
            "-q", query_genome_fp,
            "--refList", reflist_fp,
            "-o", ani_fp,
        ])

        with open(ani_fp) as f:
            for ani_result in parse_ani(f):
                yield ani_result

    def compute_ani(self, search_results):
        query = search_results[0].query
        query_fp = self.db.download_genome(query)

        subjects = list(set(r.subject for r in search_results))
        subject_fps = {self.db.download_genome(s): s for s in subjects}

        ani_results = self.get_ani(query_fp, subject_fps.keys())

        subject_results = {s: None for s in subjects}
        for res in ani_results:
            subject_fp = res["ref_fp"]
            subject = subject_fps[subject_fp]
            subject_results[subject] = res

        return [subject_results[r.subject] for r in search_results]
