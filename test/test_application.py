import collections
import os

from stackebrandtcurves.refseq import RefSeq
from stackebrandtcurves.application import StackebrandtApp

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

refseq = RefSeq(DATA_DIR)
refseq.load()

def test_main_search():
    app = StackebrandtApp(refseq)
    app.min_pctid = 95.0
    hits = app.search('GCF_001688845.2')

    pctids = {hit["sseqid"]: round(hit["pident"], 3) for hit in hits}
    assert pctids == EXPECTED_PCTIDS

def test_main_ani():
    app = StackebrandtApp(refseq)
    app.min_pctid = 95.0
    results = app.calculate_ani('GCF_001688845.2', ACCESSIONS)

    observed_anis = {a: res["ani"] for a, res in zip(ACCESSIONS, results)}
    for accession, observed_ani in observed_anis.items():
        expected_ani = EXPECTED_ANIS[accession]
        assert abs(observed_ani - expected_ani) < 0.5

def test_main_run():
    app = StackebrandtApp(refseq)
    app.min_pctid = 95.0
    results = app.run("GCF_001688845.2")
    results = list(results)

    accessions = {r.subject_seqid: r.subject_accession for r in results}
    assert accessions == EXPECTED_ACCESSIONS

    pctids = {r.subject_seqid: round(r.hit["pident"], 3) for r in results}
    assert pctids == EXPECTED_PCTIDS

    anis = {r.subject_accession: r.ani_result["ani"] for r in results}
    for accession, observed_ani in anis.items():
        expected_ani = EXPECTED_ANIS[accession]
        assert abs(observed_ani - expected_ani) < 0.5

def test_regular_search():
    app = StackebrandtApp(refseq)
    app.min_pctid = 95.0
    hits = app.regular_search("GCF_001688845.2")

    observed_pctids = {hit['sseqid']: round(hit['pident'], 3) for hit in hits}
    assert observed_pctids == EXPECTED_PCTIDS

def test_exhaustive_search():
    app = StackebrandtApp(refseq)
    app.min_pctid = 95.0
    app.max_hits = 7
    hits = app.exhaustive_search("GCF_001688845.2")

    observed_pctids = {hit['sseqid']: round(hit['pident'], 3) for hit in hits}
    assert observed_pctids == EXPECTED_PCTIDS


EXPECTED_ACCESSIONS = {
    'lcl|NZ_CP021421.1_rrna_43': 'GCF_002201515.1',
    'lcl|NZ_CP065316.1_rrna_60': 'GCF_016696845.1',
    'lcl|NZ_CP021421.1_rrna_23': 'GCF_002201515.1',
    'lcl|NZ_PUBW01000105.1_rrna_3': 'GCF_003024855.1',
    'lcl|NZ_SRYD01000003.1_rrna_36': 'GCF_004793655.1',
    'lcl|NZ_CP065316.1_rrna_40': 'GCF_016696845.1',
    'lcl|NZ_CP021421.1_rrna_34': 'GCF_002201515.1',
    'lcl|NZ_CP021421.1_rrna_38': 'GCF_002201515.1',
    'lcl|NZ_PUEE01000092.1_rrna_61': 'GCF_003024845.1',
    'lcl|NZ_CP065316.1_rrna_51': 'GCF_016696845.1',
    'lcl|NZ_CP065316.1_rrna_55': 'GCF_016696845.1',
}
    
EXPECTED_ANIS = {
    'GCF_002201515.1': 99.9958,
    'GCF_016696845.1': 99.9947,
    'GCF_003024855.1': 99.8498,
    'GCF_004793655.1': 97.8125,
    'GCF_003024845.1': 99.8737,
    'GCF_004570975.1': 79.2032,
    'GCF_004803695.1': 80.3492,
    'GCF_004793735.1': 77.9036,
    'GCF_003024925.1': 79.1948,
    'GCF_004803935.1': 81.6801,
    'GCF_004803915.1': 82.1355,
    'GCF_003024805.1': 81.8873,
    'GCF_005304985.1': 81.4932,
    'GCF_910575715.1': 81.5101,
}

ACCESSIONS = [
    'GCF_002201515.1', 'GCF_016696845.1', 'GCF_002201515.1', 'GCF_003024855.1',
    'GCF_004793655.1', 'GCF_016696845.1', 'GCF_002201515.1', 'GCF_002201515.1',
    'GCF_003024845.1', 'GCF_016696845.1', 'GCF_016696845.1',
]


EXPECTED_PCTIDS = {
    'lcl|NZ_CP021421.1_rrna_43': 100.0,
    'lcl|NZ_CP065316.1_rrna_60': 100.0,
    'lcl|NZ_CP021421.1_rrna_23': 99.74,
    'lcl|NZ_PUBW01000105.1_rrna_3': 99.74,
    'lcl|NZ_SRYD01000003.1_rrna_36': 99.74,
    'lcl|NZ_CP065316.1_rrna_40': 99.74,
    'lcl|NZ_CP021421.1_rrna_34': 99.676,
    'lcl|NZ_CP021421.1_rrna_38': 99.676,
    'lcl|NZ_PUEE01000092.1_rrna_61': 99.676,
    'lcl|NZ_CP065316.1_rrna_51': 99.676,
    'lcl|NZ_CP065316.1_rrna_55': 99.676,
}
