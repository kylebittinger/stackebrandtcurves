import logging
import os
import subprocess
import tempfile


class FastAni:
    fields = ["query_fp", "ref_fp", "ani", "fragments_aligned", "fragments_total"]
    field_types = [str, str, float, int, int]

    def __init__(self, work_dir=None):
        if work_dir is not None:
            if not os.path.exists(work_dir):
                os.makedirs(work_dir)
            self.work_dir = work_dir
        else:
            self._work_dir_obj = tempfile.TemporaryDirectory()
            self.work_dir = self._work_dir_obj.name
        logging.info(f"fastANI working directory: {self.work_dir}")

    def run(self, query_genome_fp, subject_genome_fps, threads=None):
        if threads is None:
            # https://docs.python.org/3/library/multiprocessing.html#multiprocessing.cpu_count
            threads = len(os.sched_getaffinity(0))

        reflist_fp = os.path.join(self.work_dir, "ref_list.txt")
        with open(reflist_fp, "w") as f:
            for subject_fp in subject_genome_fps:
                f.write(subject_fp)
                f.write("\n")
        ani_fp = os.path.join(self.work_dir, "ani.txt")

        args = [
            "fastANI",
            "--query",
            query_genome_fp,
            "--refList",
            reflist_fp,
            "--output",
            ani_fp,
            "--threads",
            str(threads),
            "--minFrag",
            "1",
        ]

        logging.info(f"Calling: {str(args)}")
        subprocess.check_call(args)
        logging.info("Finished fastANI call")

        with open(ani_fp) as f:
            for ani_result in self.parse(f):
                yield ani_result

    @classmethod
    def parse(cls, f):
        for line in f:
            toks = line.strip().split()
            vals = [fcn(tok) for tok, fcn in zip(toks, cls.field_types)]
            yield dict(zip(cls.fields, vals))
