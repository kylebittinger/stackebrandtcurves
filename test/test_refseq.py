import collections
import os

from stackebrandtcurves.refseq import RefSeq, parse_desc, is_low_quality

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

MockAssembly = collections.namedtuple(
    "Assembly", ["accession", "basename", "genome_url"])

a = MockAssembly(
    "GCF_001688845.2",
    "GCF_001688845.2_ASM168884v2",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/688/845/GCF_001688845.2_ASM168884v2",
)

def test_parse_desc():
    accession, attrs = parse_desc(
        "lcl|NZ_VSMC01000076.1_rrna_55 [locus_tag=FW767_RS13005] "
        "[db_xref=RFAM:RF00177] [product=16S ribosomal RNA] "
        "[location=complement(482..>596)] [gbkey=rRNA]")
    assert accession == "lcl|NZ_VSMC01000076.1_rrna_55"
    assert attrs == {
        "locus_tag": "FW767_RS13005",
        "db_xref": "RFAM:RF00177",
        "product": "16S ribosomal RNA",
        "location": "complement(482..>596)",
        "gbkey": "rRNA",
    }

def test_load_assemblies():
    db = RefSeq(DATA_DIR)
    db.load_assemblies()

    assert len(db.assemblies) == 21

    a = db.assemblies["GCF_001688845.2"]
    assert a.accession == "GCF_001688845.2"
    assert a.bioproject == "PRJNA224116"

def test_get_16S_seqs():
    db = RefSeq(DATA_DIR)
    db.load_assemblies()

    seqs = db.get_16S_seqs("GCF_001688845.2")
    seqs = list(seqs)
    assert len(seqs) == 4

    seqid, seq = seqs[0]
    assert seqid == "lcl|NZ_CP015402.2_rrna_41"
    assert seq.startswith("ACAACGAAGAGTTTGATCCTGGCTCAGGATGAACGCTAGCGACAGG")

def test_collect_genome():
    db = RefSeq(DATA_DIR)
    db.load_assemblies()

    genome_fp = db.collect_genome("GCF_001688845.2")
    assert os.path.exists(genome_fp)

def test_is_low_quality():
    assert is_low_quality("ACGTNNNNCGT", 4)
    assert not is_low_quality("ACGTNNNNCGT", 5)
