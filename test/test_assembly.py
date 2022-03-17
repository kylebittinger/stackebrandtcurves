import os

from stackebrandtcurves.assembly import RefseqAssembly

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
ASSEMBLY_SUMMARY_FP = os.path.join(DATA_DIR, "muribaculum_assembly_summary.txt")

def test_load():
    assemblies = RefseqAssembly.load(ASSEMBLY_SUMMARY_FP)
    a = assemblies["GCF_001688845.2"]
    assert a.accession == "GCF_001688845.2"
    assert a.ftp_path == (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/688/845/"
        "GCF_001688845.2_ASM168884v2")

