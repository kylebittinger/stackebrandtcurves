import collections

from stackebrandtcurves.ssu import Refseq16SDatabase

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
