import argparse
parser = argparse.ArgumentParser (description = """This program calculates a flexibility score for the amino acids of protein sequence or a family of proteins.
                                                It returns a parseable text output file and a graphical representation of the scores. """)
parser.add_argument('-i', '--input',
                    dest = "infile",
                    action = "store",
                    default = None,
                    help = "Input FASTA file or UniProt IDs")

parser.add_argument('-o', '--output',
                    dest = "outfile",
                    action = "store",
                    default = "flexibility_results",
                    help = "Flexibility analysis output files")

options = parser.parse_args()

options.infile
options.outfile
