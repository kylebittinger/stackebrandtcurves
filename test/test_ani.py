import os

from stackebrandtcurves.ani import FastAni

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

def genome_fp(fname):
    return os.path.join(DATA_DIR, 'genome_fasta', fname)

def test_fast_ani():
    query_fp = genome_fp('GCF_001688845.2_ASM168884v2_genomic.fna')
    subject_fps = [
        genome_fp('GCF_002201515.1_ASM220151v1_genomic.fna'),
        genome_fp('GCF_003024805.1_ASM302480v1_genomic.fna'),
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
