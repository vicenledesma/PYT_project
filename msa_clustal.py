#####################################
############ CLUSTAL MSA ############
#####################################

import sys
import os
from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline

def run_clustalo (input_FASTA, out_MSA_file, MSA_format = "fasta"):

    """ Wrapper to call CLUSTAL OMEGA using Python. Takes an input 
    file with FASTA sequences and creates an output file with the 
    MSA in the specified format (default: FASTA) """

    # write command

    clustalo_cline = ClustalOmegaCommandline(infile = input_FASTA, 
                                             outfile = out_MSA_file, 
                                             outfmt = MSA_format)

    # execute command

    clustalo_cline_exec = str(clustalo_cline)
    os.system(clustalo_cline_exec)  


def read_msa (MSA_file, MSA_format = "fasta"):

        """ Reads an MSA, separates it into records and returns
        it as an ID -> sequence dictionary """

        mult_aln = AlignIO.read(MSA_file, MSA_format)

        # inform of the results of the alignment

        sys.stderr.write("MSA created correctly with CLUSTAL OMEGA. ")
        sys.stderr.write(str(mult_aln) + "\n")

        # separate alignment records, create and return dict

        msa_records = {} 

        for record in mult_aln:
            msa_records[record.id] = record.seq

        return(msa_records)
