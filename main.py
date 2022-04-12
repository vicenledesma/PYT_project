
  
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
import af
import normalized_b_values
import msa_clustal
import calculate_flexibility
import output_text_graph

### Check input

if ((input_sequence and not input_family) or (not input_sequence and input_family)): # xor, only one input
    
    sys.stderr.write("Processing input...\n")

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

    try:
        input_file.write(str(consensus_seq.get_consensus_seq(input_family)))
        input_file.close()
        input_sequence = input_file.name  

    except:
        input_file.close()
        raise SystemExit("Cannot obtain consensus sequence for protein family. Are you sure you introduced a MSA in CLUSTAL format?")


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
                                            out_pssm_filename = "psiblast_uniprot_1.pssm",
                                            out_hits_filename = "psiblast_uniprot_1.out",
                                            output_format = 6) # tabular
                                    
### Get hits from PDB

AF_scores_to_download = af.select_hits_uniprot("psiblast_uniprot_1.out",
                                                number_hits = 5)

list_id_AF = af.download_hits_alphafold (AF_scores_to_download.keys(), os.getcwd() + "/PDB_downloads/")

## Extract protein sequences from the PDB files
# File for CLUSTAL MSA

hits_for_MSA_file = open('top_hits_AF.fa', 'w')

all_seq_Bnorm = {}

for AF_file in (os.listdir(os.getcwd() + "/PDB_downloads/")):
    seq_B = normalized_b_values.get_seq_B_from_PDB(os.getcwd() + "/PDB_downloads/" + AF_file)

    # Create input for MSA: templates

    hits_for_MSA_file.write (">" + AF_file.split(".")[0] + "\n")
    hits_for_MSA_file.write (seq_B["FASTA_seq"] + "\n")

    # Normalize and append to big dictionary

    seq_B_norm = normalized_b_values.normalized_b_values(seq_B)
    all_seq_Bnorm[AF_file.split(".")[0]] = seq_B_norm

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

msa_clustal.run_clustalo(hits_for_MSA_file.name, "AF_MSA_clustalo.out.fasta")

### Put indices to MSA records
## Read MSA

msa_records = msa_clustal.read_msa("AF_MSA_clustalo.out.fasta")

## Assign indices

msa_records_with_indices = msa_clustal.assign_msa_record_indices(msa_records)

## Calculate flexibility per position

sys.stderr.write("Calculating flexibility...\n")

incomplete_flexibility = calculate_flexibility.calculate_flex_pos(all_seq_Bnorm, msa_records_with_indices, AF_scores_to_download)
complete_flexibility = calculate_flexibility.complete_missing_scores(incomplete_flexibility)

### Create perseable text output file

try:
    output_text_graph.create_output_text(output_prefix, complete_flexibility)

except:
    pass

else:
    sys.stderr.write("The text output file has been created successfully!\n")

### Draw flexibility graph

try:
    output_text_graph.draw_flex_line(output_prefix, complete_flexibility)

except:
    pass

else:
    sys.stderr.write("The plot has been created successfully!\n")
