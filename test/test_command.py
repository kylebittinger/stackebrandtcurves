import collections
import os

from stackebrandtcurves.command import main, limit_results

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
    observed_output_fp = tmp_path / "output.txt"
    args = [
        "GCF_001688845.2",
        "--output-file", str(observed_output_fp),
        "--data-dir", DATA_DIR,
    ]
    main(args)
    with open(EXPECTED_OUTPUT_FP) as f:
        expected_output = f.read()
    with open(observed_output_fp) as f:
        observed_output = f.read()
    assert observed_output == expected_output
