import os

from stackebrandtcurves.refseq import RefSeq
from stackebrandtcurves.ani import FastAni

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

db = RefSeq(DATA_DIR)
db.load()

def test_fast_ani():
    query_fp = db.collect_genome('GCF_001688845.2')
    subject_fps = [
        db.collect_genome('GCF_002201515.1'),
        db.collect_genome('GCF_003024805.1'),
        ]

    app = FastAni()
    results = app.run(query_fp, subject_fps)

    for result in results:
        observed_ani = result["ani"]
        fname = os.path.basename(result["ref_fp"])
        expected_ani = EXPECTED_ANIS[fname]
        assert abs(observed_ani - expected_ani) < 0.5
    

EXPECTED_ANIS = {
    'GCF_002201515.1_ASM220151v1_genomic.fna': 99.996,
    'GCF_003024805.1_ASM302480v1_genomic.fna': 81.9096,
}
