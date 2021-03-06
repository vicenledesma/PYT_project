####################################################
############## GET HITS FROM ALPHAFOLD #############
####################################################
# LIbraries
import os
import re
import sys
from urllib import request

uniprot_file="psiblast_uniprot_1.out"

def select_hits_uniprot(blastp_outFile, number_hits=5):

    """Parses  BLASTP output file and retrieves the UniProt IDs for the hits. It selects
    only the indicated number of first hits"""

    n = 0
    UniProt_IDs = {}

    with open(blastp_outFile) as bpfile:
        for line in bpfile:
            if "100.000" not in line: # avoid repeated query protein
                m = line.split()

                if m and n<number_hits: # avoids None type objects when there are not matches
                # selects the indicated number of first hits
                    UniProt_IDs[m[1]] = m[2] # append the word that matched
                    n+=1

        if UniProt_IDs: # check that UniProt_IDs is not empty
            return(UniProt_IDs)

        else:
            raise SystemExit('No homologues found in UniProt. Exiting the program.')

def download_hits_alphafold(dict_hits_uniprot, pdb_location):
    """It downloads the PDB files of the first hits from the psi-blast from
    alphafold and saves them into the folder 'pdb'"""
    n=0

    # Create tupple list
    list_id_chain_separated = [] # initialize list for (id, chain) tupples

    if os.path.isdir(pdb_location): # checks if directory exists
        pass
    else:
        os.mkdir(pdb_location) # create directory if it does not exist

    for hit_top in dict_hits_uniprot:

    # Download PDB files

        PDB_file_location = pdb_location + hit_top + ".pdb" # directory for downloaded files

        try: # search in PDB
            request.urlretrieve('https://alphafold.ebi.ac.uk/files/AF-' + hit_top + '-F1-model_v2.pdb',PDB_file_location) # download file
            sys.stderr.write("%s PDB file downloaded successfully\n" %hit_top)
            n+=1 # number of hits found in alphafold

        except:
            pass

    print ("%s homologues found in AlphaFold" %n)
    return (list_id_chain_separated)
