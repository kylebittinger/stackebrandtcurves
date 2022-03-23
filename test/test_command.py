import os

from stackebrandtcurves.ssu import Refseq16SDatabase
from stackebrandtcurves.assembly import RefseqAssembly
from stackebrandtcurves.command import main

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
EXPECTED_OUTPUT_FP = os.path.join(DATA_DIR, "muribaculum_output.txt")

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

