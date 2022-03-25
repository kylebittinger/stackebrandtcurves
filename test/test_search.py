import collections
import os

from stackebrandtcurves.refseq import RefSeq
from stackebrandtcurves.ssu import Vsearch, count_matches, limit_hits
from stackebrandtcurves.application import StackebrandtApp

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

refseq = RefSeq(DATA_DIR)
refseq.load_assemblies()
refseq.load_seqs()

def test_count_matches():
    s1 = "ACGTNNY-G"
    s2 = "ACGGNAA-G"
    #     ^^^ ^^ ^^
    assert count_matches(s1, s2) == 7

def test_search_seq():
    app = StackebrandtApp(refseq)
    app.min_pctid = 95.0
    hits = app.regular_search("GCF_001688845.2")
    observed_pctids = {hit['sseqid']: round(hit['pident'], 2) for hit in hits}
    assert observed_pctids == EXPECTED_PCTIDS

def test_exhaustive_search():
    app = StackebrandtApp(refseq)
    app.min_pctid = 95.0
    app.max_hits = 7
    hits = app.exhaustive_search("GCF_001688845.2")
    observed_pctids = {hit['sseqid']: round(hit['pident'], 2) for hit in hits}
    assert observed_pctids == EXPECTED_PCTIDS

def test_limit_results():
    hits = [{'pident': x} for x in [90.1, 90.1, 90.1, 90.0]]
    observed = list(limit_hits(hits, 2))
    expected = [{'pident': x} for x in [90.1, 90.1, 90.0]]
    assert observed == expected

EXPECTED_PCTIDS = {
    'lcl|NZ_CP021421.1_rrna_43': 100.0,
    'lcl|NZ_CP065316.1_rrna_60': 100.0,
    'lcl|NZ_CP021421.1_rrna_23': 99.74,
    'lcl|NZ_PUBW01000105.1_rrna_3': 99.74,
    'lcl|NZ_SRYD01000003.1_rrna_36': 99.74,
    'lcl|NZ_CP065316.1_rrna_40': 99.74,
    'lcl|NZ_CP021421.1_rrna_34': 99.68,
    'lcl|NZ_CP021421.1_rrna_38': 99.68,
    'lcl|NZ_PUEE01000092.1_rrna_61': 99.68,
    'lcl|NZ_CP065316.1_rrna_51': 99.68,
    'lcl|NZ_CP065316.1_rrna_55': 99.68,
}
