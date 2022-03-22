import os

from stackebrandtcurves.assembly import RefseqAssembly

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
ASSEMBLY_SUMMARY_FP = os.path.join(DATA_DIR, "muribaculum_assembly_summary.txt")

RefseqAssembly.data_dir = DATA_DIR

def test_load_parse_summary():
    with open(ASSEMBLY_SUMMARY_FP) as f:
        assemblies = RefseqAssembly.load(f)
    a = assemblies["GCF_001688845.2"]
    assert a.accession == "GCF_001688845.2"
    assert a.ftp_path == (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/688/845/"
        "GCF_001688845.2_ASM168884v2")
    assert a.taxid == "1796646"

def test_ssu_seqs():
    a = RefseqAssembly(
        assembly_accession = "GCF_001688845.2",
        ftp_path = (
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/688/845/"
            "GCF_001688845.2_ASM168884v2"),
    )

    seqs = a.ssu_seqs
    assert len(seqs) == 4

    s0_desc, s0_seq = seqs[0]
    assert s0_desc.startswith("lcl|NZ_CP015402.2_rrna_41 ")
    assert s0_seq.startswith("ACAACGAAGAGTTTGATCCTGGCTCAGGATGAACGCTAGCGACAGG")

def test_download_genome():
    a = RefseqAssembly(
        assembly_accession = "GCF_001688845.2",
        ftp_path = (
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/688/845/"
            "GCF_001688845.2_ASM168884v2"),
    )

    a.download_genome()
    assert os.path.exists(a.genome_fp)
