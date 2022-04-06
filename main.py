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

### Check input

if ((input_sequence and not input_family) or (not input_sequence and input_family)): # xor, only one input
    sys.stderr.write("Correct input. Processing...\n")

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

### PSI-BLAST
## Get PSSM

psi_BLAST_PSSM.run_psiblast_homologues_PSSM(FASTA_seq = input_sequence, 
                                            input_iters = 5, 
                                            database = "./UniProt/UniProt.db.fasta",
                                            in_pssm_filename = None,
                                            out_pssm_filename = "psiblast_uniprot_5.pssm",
                                            out_hits_filename = "psiblast_uniprot_5.out",
                                            output_format = 6) # tabular

## Search in PDB with PSSM

psi_BLAST_PSSM.run_psiblast_homologues_PSSM(FASTA_seq = input_sequence, 
                                            input_iters = 1, 
                                            database = "./PDB_FASTA/PDB.db.fasta",
                                            in_pssm_filename = "./psiblast_uniprot_5.pssm",
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

def calculate_flex_pos(B_val_dict, pos_dict, score_dict):
    ''' Function that takes a dictionary with normalized B-values, a dictionary 
    with the indexed position of a MSA and a dictionary with BLAST scores and returns
    a list of (residue, flexibility) tuples for a query. '''
    print (B_val_dict)
    flexibility_position_list = []

    query_position = pos_dict.pop("query")
    #print(query_position)

    for query_index_tuple in query_position:

        B_vals_weighted = []
        scores = []
        
        query_residue = query_index_tuple[0]
        query_MSA_index = query_index_tuple[1]

        for template_ID, template_index_list in pos_dict.items():
            for template_index_tuple in template_index_list:

                template_MSA_index = template_index_tuple[1]
                template_PDB_index = template_index_tuple[2]

                if query_MSA_index == template_MSA_index:
                    B_val_raw = float(B_val_dict[template_ID]["B_val_list"][template_PDB_index][1])
                    B_vals_weighted.append(B_val_raw * float(score_dict[template_ID]))
                    scores.append(float(score_dict[template_ID]))
                    print(B_vals_weighted)
                    print(scores)
                       
        if (B_vals_weighted and scores):
            flexibility_position_list.append((query_residue,(sum(B_vals_weighted)/sum(scores))))
                
    print(flexibility_position_list)

    
        




### Pruebas

calculate_flex_pos(all_seq_Bnorm, msa_records_with_indices, PDB_scores_to_download)


