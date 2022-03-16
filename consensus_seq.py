import sys
from Bio import AlignIO
from Bio.Align import AlignInfo
def consensus_seq(fasta_file):
    """Creates a consensus sequence
    from an input MSA file in clustal format.
    If the percentage of the most common residue type is
    greater than 70%, then we will add that
    residue type; otherwise an X will be added."""
    alignment = AlignIO.read(fasta_file, 'clustal')
    summary_align = AlignInfo.SummaryInfo(alignment)
    con_seq = summary_align.dumb_consensus()
    return (con_seq)
fam_cons_seq = consensus_seq('alignment_glo.aln')
print(fam_cons_seq)
