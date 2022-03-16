#########################################
############ PSI-BLAST + IDs ############
#########################################

import os
import re
from Bio.Blast.Applications import NcbipsiblastCommandline
from Bio.Blast import Record
from Bio.Blast import NCBIXML
from Bio import SearchIO


def run_psiblast_homologues_PSSM (FASTA_seq, 
                                  input_iters,
                                  database,
                                  in_pssm_filename,
                                  out_pssm_filename,
                                  out_hits_filename,
                                  output_format = 6
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
                                                 out = out_hits_filename,
                                                 outfmt = output_format)
                                            
    else: # no input PSSM provided
        psiblast_cline = NcbipsiblastCommandline(query = FASTA_seq, 
                                                 num_iterations = input_iters,
                                                 db = database,
                                                 out_pssm = out_pssm_filename,
                                                 out = out_hits_filename,
                                                 outfmt = output_format)

    # execute command
    
    psiblast_cline_exec = str(psiblast_cline)
    os.system(psiblast_cline_exec)        


def get_IDs_from_blastp_PDB (blastp_outFile):

    """ Parses  BLASTP output file and retrieves the PDB IDs and the scores
    for the hits. Returns them as a dictionary: ID -> score """

    # regular expresion for PDB IDs

    PDB_ID_score = {}

    with open(blastp_outFile) as bpfile:
        for line in bpfile:

            records = line.split()
            PDB_ID = records[1]
            score = records[11]

            PDB_ID_score[PDB_ID] = score

        if PDB_ID_score: # check that UniProt_IDs is not empty
            return(PDB_ID_score)  

        else:
            raise SystemExit('No homologues found in UniProt. Exiting the program.')    

