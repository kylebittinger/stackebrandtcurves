import collections
import os

from stackebrandtcurves.refseq import RefSeq

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
ASSEMBLY_SUMMARY_FP = os.path.join(DATA_DIR, "muribaculum_assembly_summary.txt")

MockAssembly = collections.namedtuple(
    "Assembly", ["accession", "basename", "genome_url"])

a = MockAssembly(
    "GCF_001688845.2",
    "GCF_001688845.2_ASM168884v2",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/688/845/GCF_001688845.2_ASM168884v2",
)

def test_ssu_seqs():
    db = RefSeq(DATA_DIR)
    rna_fp = db.download_rna(a)
    assert rna_fp == os.path.join(
        DATA_DIR, "rna_fasta", "GCF_001688845.2_ASM168884v2_rna_from_genomic.fna")

    seqs = db.get_16S_seqs(a)
    seqs = list(seqs)
    assert len(seqs) == 4

    desc, seq = seqs[0]
    assert desc.startswith("lcl|NZ_CP015402.2_rrna_41 ")
    assert seq.startswith("ACAACGAAGAGTTTGATCCTGGCTCAGGATGAACGCTAGCGACAGG")

def test_download_genome():
    db = RefSeq(DATA_DIR)
    genome_fp = db.download_genome(a)
    assert os.path.exists(genome_fp)
