import os
import subprocess
import tempfile

from .parse import parse_ani

class FastAni:
    def __init__(self, work_dir=None):
        if work_dir is not None:
            if not os.path.exists(work_dir):
                os.makedirs(work_dir)
            self.work_dir = work_dir
        else:
            self._work_dir_obj = tempfile.TemporaryDirectory()
            self.work_dir = self._work_dir_obj.name

    def run(self, query_genome_fp, subject_genome_fps):
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
