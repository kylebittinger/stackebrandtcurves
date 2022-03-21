import os

from stackebrandtcurves.ssu import Refseq16SDatabase
from stackebrandtcurves.assembly import RefseqAssembly
from stackebrandtcurves.command import main

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
ASSEMBLY_SUMMARY_FP = os.path.join(DATA_DIR, "muribaculum_assembly_summary.txt")
TEST_FASTA = os.path.join(DATA_DIR, "muribaculum.fasta")
TEST_ACCESSIONS = os.path.join(DATA_DIR, "muribaculum_accessions.txt")
OUTPUT_FP = os.path.join(DATA_DIR, "muribaculum_output.txt")

RefseqAssembly.summary_fp = ASSEMBLY_SUMMARY_FP
RefseqAssembly.rna_dir = DATA_DIR
RefseqAssembly.genome_dir = DATA_DIR

def test_main(tmpdir):
    observed_output_fp = os.path.join(DATA_DIR, "observed_test_output.txt")
    main(["GCF_001688845.2", "--output-file", observed_output_fp])
    with open(OUTPUT_FP) as f:
        expected_output = f.read()
    with open(observed_output_fp) as f:
        observed_output = f.read()
    assert observed_output == expected_output

