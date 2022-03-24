from io import StringIO

OUTPUT_FIELDS = [
    "query_assembly", "subject_assembly", "query_seqid", "subject_seqid",
    "pctid", "ani", "fragments_aligned", "fragments_total"
    ]
OUTPUT_TYPES = [str, str, str, str, float, float, int, int]

def parse_output(f):
    next(f) # skip header
    for line in f:
        toks = line.strip().split("\t")
        vals = [fcn(tok) for tok, fcn in zip(toks, OUTPUT_TYPES)]
        yield dict(zip(OUTPUT_FIELDS, vals))

ANI_FIELDS = [
    "query_fp", "ref_fp", "ani", "fragments_aligned", "fragments_total"]
ANI_TYPES = [str, str, float, int, int]

def parse_pairwise_ani(f):
    for line in f:
        toks = line.strip().split()
        vals = [fcn(tok) for tok, fcn in zip(toks, ANI_TYPES)]
        yield dict(zip(ANI_FIELDS, vals))

ASSEMBLY_SUMMARY_COLS = [
    "assembly_accession", "bioproject", "biosample", "wgs_master",
    "refseq_category", "taxid", "species_taxid", "organism_name",
    "infraspecific_name", "isolate", "version_status",
    "assembly_level", "release_type", "genome_rep", "seq_rel_date",
    "asm_name", "submitter", "gbrs_paired_asm", "paired_asm_comp",
    "ftp_path", "excluded_from_refseq", "relation_to_type_material"
]

def parse_assembly_summary(f):
    for line in f:
        line = line.rstrip("\n")
        if line.startswith("#") or (line == ""):
            continue
        toks = line.split("\t")
        vals = dict(zip(ASSEMBLY_SUMMARY_COLS, toks))
        if vals["ftp_path"] == "na":
            continue
        yield vals

def parse_species_names(f):
    for desc, seq in parse_fasta(f):
        vals = desc.split("\t", maxsplit=1)
        accession = vals[0]
        if len(vals) == 2:
            species_name = vals[1]
        else:
            species_name = accession
        yield accession, species_name

def parse_fasta(f, trim_desc=False):
    """Parse a FASTA format file.

    Parameters
    ----------
    f : File object or iterator returning lines in FASTA format.

    Returns
    -------
    An iterator of tuples containing two strings
        First string is the sequence description, second is the
        sequence.

    Notes
    -----
    This function removes whitespace in the sequence and translates
    "U" to "T", in order to accommodate FASTA files downloaded from
    SILVA and the Living Tree Project.
    """
    f = iter(f)
    try:
        desc = next(f).strip()[1:]
        if trim_desc:
            desc = desc.split()[0]
    except StopIteration:
        return
    seq = StringIO()
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            yield desc, seq.getvalue()
            desc = line[1:]
            if trim_desc:
                desc = desc.split()[0]
            seq = StringIO()
        else:
            seq.write(line.replace(" ", "").replace("U", "T"))
    yield desc, seq.getvalue()


def write_fasta(f, seqs):
    for desc, seq in seqs:
        f.write(">{0}\n{1}\n".format(desc, seq))


def load_fasta(filepath, trim_desc=True):
    """Load all sequences from a FASTA file

    Parameters
    ----------
    fasta_fp : Input filepath, FASTA format.

    Returns
    -------
    A dictionary mapping sequence identifiers to sequences.
    """
    with open(filepath) as f:
        seqs = parse_fasta(f, trim_desc=trim_desc)
        return dict(seqs)


def parse_greengenes_accessions(f):
    for line in f:
        if line.startswith("#"):
            continue
        line = line.strip()
        yield line.split("\t")

def parse_results(f):
    float_fields = ["probability_incompatible"]
    int_fields = ["region_mismatches", "region_positions"]
    header_line = next(f)
    header_line = header_line.rstrip()
    fields = header_line.split("\t")
    for line in f:
        line = line.rstrip()
        vals = line.split("\t")
        res = dict(zip(fields, vals))
        for field, val in res.items():
            if field in float_fields:
                res[field] = float(val)
            elif field in int_fields:
                res[field] = int(val)
        yield res
