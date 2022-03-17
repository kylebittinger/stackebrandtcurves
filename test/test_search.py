import collections
import os

from stackebrandtcurves.ssu import Refseq16SDatabase
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

MockAssembly = collections.namedtuple("Assembly", ["accession", "ssu_seqs"])

def test_add_assembly():
    a_seqs = [
        ("seq1", "GCTCGCATCGAT"),
        ("seq2", "GCTCGCATCGAT"), # duplicate, will not be added
        ("seq3", "TGCTCAGTCGT"),
    ]
    a = MockAssembly("GCF_001688845.2", a_seqs)
    db = Refseq16SDatabase()
    db.add_assembly(a)

    assert db.assemblies["seq1"] == a
    assert "seq2" not in db.assemblies
    assert db.assemblies["seq3"] == a

    assert db.seqs["seq1"] == "GCTCGCATCGAT"
    assert "seq2" not in db.seqs
    assert db.seqs["seq3"] == "TGCTCAGTCGT"

    assert db.seqids_by_assembly["GCF_001688845.2"] == ["seq1", "seq3"]
