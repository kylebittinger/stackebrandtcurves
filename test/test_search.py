import collections
import os

from stackebrandtcurves.refseq import RefSeq
from stackebrandtcurves.ssu import count_matches, limit_hits

def test_count_matches():
    s1 = "ACGTNNY-G"
    s2 = "ACGGNAA-G"
    #     ^^^ ^^ ^^
    assert count_matches(s1, s2) == 7

def test_limit_hits():
    hits = [{'pident': x} for x in [90.1, 90.1, 90.1, 90.0]]
    observed = list(limit_hits(hits, 2))
    expected = [{'pident': x} for x in [90.1, 90.1, 90.0]]
    assert observed == expected
