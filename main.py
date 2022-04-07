#########################################
############## MAIN PROGRAM #############
#########################################
### Import modules

from ast import dump
import sys
import os
from argparser import *
import consensus_seq
import psi_BLAST_PSSM
import top_hits_pdb
import normalized_b_values
import msa_clustal
import calculate_flexibility

### Check input

if ((input_sequence and not input_family) or (not input_sequence and input_family)): # xor, only one input
    
    sys.stderr.write("Correct input. Processing...\n")

    # remove directory if it already exists

    if os.path.isdir("flexibility"):
        os.system("rm -r flexibility")

    # create new directory for results
    
    os.mkdir("flexibility")

    # move input to new directory

    if input_sequence:
        os.system("cp " + input_sequence + " ./flexibility")

    else: 
         os.system("cp " + input_family + " ./flexibility")

    os.chdir("./flexibility")


elif (input_sequence and input_family):
    raise SystemExit("Incorrect input. Please, make sure you introduce only one type of input (protein sequence or MSA of protein family in CLUSTAL format).")

else:
    raise SystemExit("No input. Please, introduce an input option.")

### The input is a MSA of a protein family

if input_family: # convert protein family to consensus sequence
    input_file = open('consensus_seq.fa', 'w')
    input_file.write(">query\n")
    input_file.write(str(consensus_seq.get_consensus_seq(input_family)))
    input_file.close()
    input_sequence = input_file.name

### PSI-BLAST or BLASTP
# ## Get PSSM

# psi_BLAST_PSSM.run_psiblast_homologues_PSSM(FASTA_seq = input_sequence, 
#                                             input_iters = 1, 
#                                             database = "../UniProt/UniProt.db.fasta",
#                                             in_pssm_filename = None,
#                                             out_pssm_filename = "psiblast_uniprot_5.pssm",
#                                             out_hits_filename = "psiblast_uniprot_5.out",
#                                             output_format = 6) # tabular

## Search in UniProt with PSSM

psi_BLAST_PSSM.run_psiblast_homologues_PSSM(FASTA_seq = input_sequence, 
                                            input_iters = 1, 
                                            database = "../UniProt/UniProt.db.fasta",
                                            in_pssm_filename = None,
                                            out_pssm_filename = "psiblast_pdb_1.pssm",
                                            out_hits_filename = "psiblast_pdb_1.out",
                                            output_format = 6) # tabular

### Get hits from PDB

PDB_scores_to_download = top_hits_pdb.get_IDs_from_blastp_PDB("psiblast_pdb_1.out",
                                                        number_hits = 7)

list_id_PDB = top_hits_pdb.download_pdb (PDB_scores_to_download.keys(), os.getcwd() + "/PDB_downloads/")

top_hits_pdb.select_chain_from_pdb(list_id_PDB, os.getcwd() + "/PDB_downloads/", os.getcwd() + "/PDB_downloads/split/")

## Extract protein sequences from the PDB files
# File for CLUSTAL MSA

hits_for_MSA_file = open('top_hits_split_chain.fa', 'w')

all_seq_Bnorm = {}

for PDB_split_file in (os.listdir(os.getcwd() + "/PDB_downloads/split")):
    seq_B = normalized_b_values.get_seq_B_from_PDB(os.getcwd() + "/PDB_downloads/split/" + PDB_split_file)

    # Create input for MSA: templates

    hits_for_MSA_file.write (">" + PDB_split_file.split(".")[0] + "\n")
    hits_for_MSA_file.write (seq_B["FASTA_seq"] + "\n")

    # Normalize and append to big dictionary

    seq_B_norm = normalized_b_values.normalized_b_values(seq_B)
    all_seq_Bnorm[PDB_split_file.split(".")[0]] = seq_B_norm

# Add query to MSA input

hits_for_MSA_file.write(">query\n")

if input_family:

    hits_for_MSA_file.write(str(consensus_seq.get_consensus_seq(input_family)))

elif input_sequence:
    fd = open (input_sequence, "r")
    for line in fd:
        if not line.startswith(">"):
            hits_for_MSA_file.write(line)

    fd.close()

hits_for_MSA_file.close()

### CLUSTAL MSA

msa_clustal.run_clustalo(hits_for_MSA_file.name, "PDB_MSA_clustalo.out.fasta")

### Put indices to MSA records
## Read MSA

msa_records = msa_clustal.read_msa("PDB_MSA_clustalo.out.fasta")

## Assign indices

msa_records_with_indices = msa_clustal.assign_msa_record_indices(msa_records)

## Calculate flexibility per position

incomplete_flexibility = calculate_flexibility.calculate_flex_pos(all_seq_Bnorm, msa_records_with_indices, PDB_scores_to_download)
complete_flexibility = calculate_flexibility.complete_missing_scores(incomplete_flexibility)

### Create perseable text output file

def create_output_file (prefix, complete_flex_list):

    '''Creates parseable text file with the flexibility scores.'''

    fd = open(prefix + ".txt", "w")
    fd.write("Residue\t\t\tFlexibility\t\tConfidence\n")

    for residue in complete_flex_list:
        fd.write(residue[0] + "\t\t\t" + "%+.3f\t\t\t" % (residue[1]) + str(residue[2]) + "\n")

    fd.close()

create_output_file(output_prefix, complete_flexibility)