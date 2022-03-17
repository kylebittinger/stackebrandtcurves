import os

from stackebrandtcurves.assembly import RefseqAssembly

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
ASSEMBLY_SUMMARY_FP = os.path.join(DATA_DIR, "muribaculum_assembly_summary.txt")

RefseqAssembly.summary_fp = ASSEMBLY_SUMMARY_FP

def test_load_parse_summary():
    assemblies = RefseqAssembly.load()
    a = assemblies["GCF_001688845.2"]
    assert a.accession == "GCF_001688845.2"
    assert a.ftp_path == (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/688/845/"
        "GCF_001688845.2_ASM168884v2")
    assert a.taxid == "1796646"

