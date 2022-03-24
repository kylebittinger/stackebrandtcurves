import collections
import os

from stackebrandtcurves.command import main, limit_results
from stackebrandtcurves.parse import parse_output

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
EXPECTED_OUTPUT_FP = os.path.join(DATA_DIR, "muribaculum_output.txt")

MockResult = collections.namedtuple("SearchResult", ["pctid"])

def test_limit_results():
    pctids = [90.1, 90.1, 90.1, 90.0]
    results = [MockResult(x) for x in pctids]
    limited_results = list(limit_results(results, 2))
    print(limited_results)
    print(results)
    assert limited_results == results[1:]


def test_main(tmp_path):
    output_fp = tmp_path / "output.txt"
    args = [
        "GCF_001688845.2",
        "--output-file", str(output_fp),
        "--data-dir", DATA_DIR,
    ]
    main(args)
    with open(output_fp) as f:
        output = list(parse_output(f))
    pctids = {x["subject_seqid"]: x["pctid"] for x in output}
    assert pctids == EXPECTED_PCTIDS

    anis = {x["subject_seqid"]: x["ani"] for x in output}
    for subj, observed_ani in anis.items():
        expected_ani = EXPECTED_ANIS[subj]
        assert abs(expected_ani - observed_ani) < 1

    subjects = set(x["subject_seqid"] for x in output)
    assert subjects == set(EXPECTED_ANIS.keys())

    assemblies = set(x["subject_assembly"] for x in output)
    assert assemblies == set(EXPECTED_ASSEMBLIES)


EXPECTED_ASSEMBLIES = [
    'GCF_003024845.1', 'GCF_910575715.1', 'GCF_004793735.1', 'GCF_016696845.1',
    'GCF_004803915.1', 'GCF_004793655.1', 'GCF_003024925.1', 'GCF_004570975.1',
    'GCF_003024805.1', 'GCF_002201515.1', 'GCF_005304985.1', 'GCF_004803935.1',
    'GCF_003024855.1', 'GCF_004803695.1',
]

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
    'lcl|NZ_SPPC01000028.1_rrna_28': 92.21,
    'lcl|NZ_CP039393.1_rrna_13': 92.02,
    'lcl|NZ_CP039393.1_rrna_61': 92.02,
    'lcl|NZ_SRYY01000027.1_rrna_14': 91.95,
    'lcl|NZ_CP039393.1_rrna_39': 91.89,
    'lcl|NZ_CP039393.1_rrna_49': 91.89,
    'lcl|NZ_PUED01000303.1_rrna_29': 90.72,
    'lcl|NZ_CP039547.1_rrna_21': 90.66,
    'lcl|NZ_CP039396.1_rrna_31': 90.6,
    'lcl|NZ_CP039396.1_rrna_13': 90.53,
    'lcl|NZ_PUEC01000003.1_rrna_37': 90.53,
    'lcl|NZ_CP039547.1_rrna_35': 90.53,
    'lcl|NZ_CP040121.1_rrna_39': 90.53,
    'lcl|NZ_CP040121.1_rrna_53': 90.53,
    'lcl|NZ_CP039396.1_rrna_26': 90.47,
    'lcl|NZ_CP039396.1_rrna_53': 90.46,
    'lcl|NZ_CP039547.1_rrna_15': 90.46,
    'lcl|NZ_CP039547.1_rrna_65': 90.46,
    'lcl|NZ_CP040121.1_rrna_11': 90.46,
    'lcl|NZ_CP040121.1_rrna_33': 90.46,
    'lcl|NZ_CAJTGC010000093.1_rrna_72': 90.33,
}

EXPECTED_ANIS = {
    'lcl|NZ_CP021421.1_rrna_43': 99.9957,
    'lcl|NZ_CP065316.1_rrna_60': 99.9948,
    'lcl|NZ_CP021421.1_rrna_23': 99.9957,
    'lcl|NZ_PUBW01000105.1_rrna_3': 99.869,
    'lcl|NZ_SRYD01000003.1_rrna_36': 97.7767,
    'lcl|NZ_CP065316.1_rrna_40': 99.9948,
    'lcl|NZ_CP021421.1_rrna_34': 99.9957,
    'lcl|NZ_CP021421.1_rrna_38': 99.9957,
    'lcl|NZ_PUEE01000092.1_rrna_61': 99.8737,
    'lcl|NZ_CP065316.1_rrna_51': 99.9948,
    'lcl|NZ_CP065316.1_rrna_55': 99.9948,
    'lcl|NZ_SPPC01000028.1_rrna_28': 79.2032,
    'lcl|NZ_CP039393.1_rrna_13': 80.3873,
    'lcl|NZ_CP039393.1_rrna_61': 80.3873,
    'lcl|NZ_SRYY01000027.1_rrna_14': 77.9219,
    'lcl|NZ_CP039393.1_rrna_39': 80.3873,
    'lcl|NZ_CP039393.1_rrna_49': 80.3873,
    'lcl|NZ_PUED01000303.1_rrna_29': 79.1306,
    'lcl|NZ_CP039547.1_rrna_21': 81.5896,
    'lcl|NZ_CP039396.1_rrna_31': 82.1546,
    'lcl|NZ_CP039396.1_rrna_13': 82.1546,
    'lcl|NZ_PUEC01000003.1_rrna_37': 81.9231,
    'lcl|NZ_CP039547.1_rrna_35': 81.5896,
    'lcl|NZ_CP040121.1_rrna_39': 81.5162,
    'lcl|NZ_CP040121.1_rrna_53': 81.5162,
    'lcl|NZ_CP039396.1_rrna_26': 82.1546,
    'lcl|NZ_CP039396.1_rrna_53': 82.1546,
    'lcl|NZ_CP039547.1_rrna_15': 81.5896,
    'lcl|NZ_CP039547.1_rrna_65': 81.5896,
    'lcl|NZ_CP040121.1_rrna_11': 81.5162,
    'lcl|NZ_CP040121.1_rrna_33': 81.5162,
    'lcl|NZ_CAJTGC010000093.1_rrna_72': 81.5595,
}
