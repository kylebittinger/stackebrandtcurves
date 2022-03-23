import collections
import os

from stackebrandtcurves.ssu import Refseq16SDatabase
from stackebrandtcurves.assembly import RefseqAssembly

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
ASSEMBLY_SUMMARY_FP = os.path.join(DATA_DIR, "muribaculum_assembly_summary.txt")
TEST_FASTA = os.path.join(DATA_DIR, "muribaculum.fasta")
TEST_ACCESSIONS = os.path.join(DATA_DIR, "muribaculum_accessions.txt")

MockAssembly = collections.namedtuple("Assembly", ["accession", "ssu_seqs"])

class MockRefSeq:
    data_dir = DATA_DIR
    def get_16S_seqs(self, assembly):
        return [
            ("seq1", "GCTCGCATCGAT"),
            ("seq2", "GCTCGCATCGAT"), # duplicate, will not be added
            ("seq3", "TGCTCAGTCGT"),
        ]

def test_search_seq():
    db = Refseq16SDatabase(MockRefSeq())
    with open(ASSEMBLY_SUMMARY_FP) as f:
        assemblies = RefseqAssembly.parse(f)
        db.load({a.accession: a for a in assemblies})
    hits = db.search_seq("lcl|NZ_CP015402.2_rrna_41", TEST_SEQ, min_pctid = 95.0)
    hits = list(hits)
    assert len(hits) == 13
    last_hit = hits[-1]
    assert last_hit.subject.accession == "GCF_004792675.1"

def test_exhaustive_search():
    db = Refseq16SDatabase(MockRefSeq())
    with open(ASSEMBLY_SUMMARY_FP) as f:
        assemblies = RefseqAssembly.parse(f)
        db.load({a.accession: a for a in assemblies})
    hits = db.exhaustive_search(
        "lcl|NZ_CP015402.2_rrna_41", TEST_SEQ, min_pctid = 95.0)
    hits = list(hits)
    assert len(hits) == 13
    last_hit = hits[-1]
    assert last_hit.subject.accession == "GCF_004792675.1"

def test_add_assembly():
    a = MockAssembly("GCF_001688845.2", None)
    db = Refseq16SDatabase(MockRefSeq())
    db.add_assembly(a)

    assert db.assemblies["seq1"] == a
    assert "seq2" not in db.assemblies
    assert db.assemblies["seq3"] == a

    assert db.seqs["seq1"] == "GCTCGCATCGAT"
    assert "seq2" not in db.seqs
    assert db.seqs["seq3"] == "TGCTCAGTCGT"

    assert db.seqids_by_assembly["GCF_001688845.2"] == ["seq1", "seq3"]


TEST_SEQ = (
    "ACAACGAAGAGTTTGATCCTGGCTCAGGATGAACGCTAGCGACAGGCCTAACACATGCAAGTCGAGGGGCAGCGGGGAGC"
    "GAAGCTTGCTTTGCTCCGCCGGCGACCGGCGCACGGGTGAGTAACGCGTATGCAACCTGCCCGTTTCAGCGGGATAACCC"
    "GGAGAAATCCGGCCTAATACCGCATGTGGCCGAGGGGGAGGCATCTTCTTTCGGCCAAAGGAGGCAACTCCGGAGACGGA"
    "TGGGCATGCGTGACATTAGCTTGTTGGCGGGGTAACGGCCCACCAAGGCGACGATGTCTAGGGGTTCTGAGAGGAAGGTC"
    "CCCCACACTGGTACTGAGACACGGACCAGACTCCTACGGGAGGCAGCAGTGAGGAATATTGGTCAATGGGCGAGAGCCTG"
    "AACCAGCCAAGTCGCGTGAAGGATGACGGCCCTACGGGTTGTAAACTTCTTTTGTCAGGGAGCAATTGAGTCCACGTGTG"
    "GGCTTAGCGAGAGTACCTGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTT"
    "ATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGAATATCAAGTCAGCGGTAAAAATTCGGGGCTCAACCCCGTCGTG"
    "CCGTTGAAACTGATGTTCTTGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCA"
    "GAACTCCGATTGCGAAGGCAGCATACCGGCGCCCAACTGACGCTGAAGCACGAAAGCGTGGGTATCGAACAGGATTAGAT"
    "ACCCTGGTAGTCCACGCAGTAAACGATGAATGCTAATTGTCCGGAGGGATTGACCTCTGGGTGATACAGCGAAAGCGTTA"
    "AGCATTCCACCTGGGGAGTACGCCGGCAACGGTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGAGGAACATGT"
    "GGTTTAATTCGATGATACGCGAGGAACCTTACCCGGGCTCAAACGCAGGAGGAATACATTTGAAAGGGTGTAGCTCTACG"
    "GAGTCTTCTGCGAGGTGCTGCATGGTTGTCGTCAGCTCGTGCCGTGAGGTGTCGGCTTAAGTGCCATAACGAGCGCAACC"
    "CCTATCGACAGTTGCCAACAGGTTAAGCTGGGAACTCTGTCGAGACTGCCGGCGCAAGCTGAGAGGAAGGCGGGGATGAC"
    "GTCAAATCAGCACGGCCCTTACGTCCGGGGCGACACACGTGTTACAATGGCGGCCACAGCGGGAAGCCAGGCGGCGACGC"
    "CGAGCGGAACCCGAAAAGCCGTCTCAGTTCGGATTGGAGTCTGCAACCCGACTCCATGAAGCTGGATTCGCTAGTAATCG"
    "CGCATCAGCCACGGCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCAAGCCATGGGAGCCGGGAGTGCCTG"
    "AAGTACGTGACCGCAAGGAGCGTCCTAGGGTAAGACCGGTGACTGGGGCTAAGTCGTAACAAGGTAGCCGTACCGGAAGG"
    "TGCGGCTGGAACACCTCCTTT"
)
