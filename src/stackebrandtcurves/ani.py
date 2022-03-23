import argparse
import collections
import itertools
import os
import os.path
import random
import re
import shutil
import subprocess
import urllib.error

from .parse import parse_pairwise_ani

class AniApplication:
    def __init__(self, db):
        self.db = db
        self.data_dir = self.db.data_dir
        self.ani_cache = {}

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
        ani_fp = os.path.join(self.data_dir, "temp_ani.txt")
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
