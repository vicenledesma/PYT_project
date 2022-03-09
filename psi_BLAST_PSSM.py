#########################################
############ PSI-BLAST + IDs ############
#########################################

import os
import re
from Bio.Blast.Applications import NcbipsiblastCommandline

def run_psiblast_homologues_PSSM (FASTA_seq, 
                                  input_iters,
                                  database,
                                  in_pssm_filename,
                                  out_pssm_filename,
                                  out_hits_filename
                                  ):

    """ Wrapper for running PSI-BLAST from the command line. It takes a sequence 
    in FASTA format and returns an output file with the results of the search in the
    database and the PSSM corresponding to the MSA with the hits """

    # write command

    if in_pssm_filename: #  input PSSM provided
        psiblast_cline = NcbipsiblastCommandline(num_iterations = input_iters,
                                                 db = database,
                                                 in_pssm = in_pssm_filename,
                                                 out_pssm = out_pssm_filename,
                                                 out = out_hits_filename)
                                            
    else: # no input PSSM provided
        psiblast_cline = NcbipsiblastCommandline(query = FASTA_seq, 
                                                 num_iterations = input_iters,
                                                 db = database,
                                                 out_pssm = out_pssm_filename,
                                                 out = out_hits_filename)

    # execute command
    
    psiblast_cline_exec = str(psiblast_cline)
    os.system(psiblast_cline_exec)        


def get_IDs_from_blastp_PDB (blastp_outFile):

    """ Parses  BLASTP output file and retrieves the PDB IDs for the hits """

    # regular expresion for PDB IDs
    PDB_regex = re.compile("[0-9][a-zA-Z_0-9]{3}_[A-Z]")

    PDB_IDs = set() # avoid repetition

    with open(blastp_outFile) as bpfile:
        for line in bpfile:
            m = PDB_regex.match(line) # match

            if m: # avoids None type objects when there are not matches
                PDB_IDs.add(m.group()) # append the word that matched
        
        if PDB_IDs: # check that UniProt_IDs is not empty
            print (PDB_IDs)
            return(PDB_IDs)   

        else:
            raise SystemExit('No homologues found in UniProt. Exiting the program.')    
