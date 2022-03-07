#########################################
############ PSI-BLAST + IDs ############
#########################################

import os
from Bio.Blast.Applications import NcbipsiblastCommandline
import re

def run_psiblast_homologues(FASTA_seq):

    """ Wrapper for running PSI-BLAST from the command line. It takes a sequence 
    in FASTA format and returns an output file with the results of the search in 
    Uniprot (with 5 iterations) """

    # write command

    psiblast_cline = NcbipsiblastCommandline(query = FASTA_seq, 
                                             num_iterations = 5,
                                             db = "/home/vicente/Desktop/UniProt/UniProt.db.fasta",
                                             # out_pssm = "query_uniprot_matrix.pssm",
                                             out = "query_uniprot_hits.out")
    # execute command

    psiblast_cline_exec = str(psiblast_cline)
    os.system(psiblast_cline_exec)        


def get_IDs_from_psiblast (psiblast_outFile):

    """ Parses PSI_BLAST output file and retrieves the UniProt IDs for
    the hits in the last (fifth) iteration """

    # regular expresion for UniProt IDs
    UP_regex = re.compile("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")

    UniProt_IDs = set() # avoid repetition
    start_reading = False

    with open(psiblast_outFile) as psifile:
        for line in psifile:
            if line.strip() == 'Results from round 5': # start reading in the last iteration
                start_reading = True

            elif start_reading:
                m = UP_regex.match(line) # match
                if m: # avoids None type objects when there are not matches
                    UniProt_IDs.add(m.group()) # append the word that matched
        
        if UniProt_IDs: # check that UniProt_IDs is not empty
            return(UniProt_IDs)   

        else:
            raise SystemExit('No homologues found in UniProt. Exiting the program.')    


