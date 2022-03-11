import os
import re
import sys
from urllib import request
from Bio import PDB

pdb_file="outPDB.hits"
def get_IDs_from_blastp_PDB (blastp_outFile, number_hits=20):
    """ Parses  BLASTP output file and retrieves the PDB IDs for the hits. It selects
    only the indicated number of first hits """
    n=0

    # regular expresion for PDB IDs
    PDB_regex = re.compile("[0-9][a-zA-Z_0-9]{3}_[A-Z]")

    PDB_IDs = set() # avoid repetition

    with open(blastp_outFile) as bpfile:
        for line in bpfile:
            m = PDB_regex.match(line) # match

            if m and n<number_hits: # avoids None type objects when there are not matches
            # selects the indicated number of first hits
                PDB_IDs.add(m.group()) # append the word that matched
                n=n+1
        if PDB_IDs: # check that UniProt_IDs is not empty
            #print (PDB_IDs)
            return(PDB_IDs)

        else:
            raise SystemExit('No homologues found in UniProt. Exiting the program.')

list_id_chain=get_IDs_from_blastp_PDB(pdb_file,2)

dir_pdb=("/Users/sarapoloalonso/Documents/Apuntes_2/PYT/Proyecto/pdb_down/")


def download_pdb (list_id_chain, pdb_location):
    """Separates the PDB ids and chains, then downloads the PDB files
    of the 20 first hits of the psi-blast and finally saves them in the folder 'pdb' """
    list_id_chain_separated=[] # initialize
    if os.path.isdir(pdb_location): # checks if directory exists
        pass
    else:
        os.mkdir(pdb_location) # create directory if it does not exist

    for hit_20 in list_id_chain:
        PDB_all=hit_20.split('_') # separate id and chains
        id=PDB_all[0].upper() # pdb id
        chain=PDB_all[1].upper() # pdb homolog chain
        tuple_id_chain=(id, chain)
        list_id_chain_separated.append(tuple_id_chain) #list with tupples
        # id as first element, chain as second element
        PDB_file_location=pdb_location+id+'.pdb' # directory for downloaded files
        try:
            # Download pdb files
            request.urlretrieve('https://files.rcsb.org/download/'+id+'.pdb', PDB_file_location) #download file
            sys.stderr.write("%s PDB file downloaded successfully\n" %id)
        except:
            raise SystemExit("ID PDB %s was not found." %id)
    return (list_id_chain_separated)

list_chain_separated=download_pdb(list_id_chain, dir_pdb)

dir_split='/Users/sarapoloalonso/Documents/Apuntes_2/PYT/Proyecto/pdb_down/split/'
def select_chain_from_pdb(id_chain_pdb, pdb_location, pdb_split_location):
    """Split the PDB files, creating new ones that
    contain only the data of the selected chains"""
    if os.path.isdir(pdb_split_location): # checks if directory exists
        pass
    else:
        os.mkdir(pdb_split_location) # creates it if it does not exist

    for file_name in id_chain_pdb:
        current_file=pdb_location+file_name[0]+'.pdb' # name of downloaded file
        selected_chain=file_name[1] # chain
        output_file=pdb_split_location+file_name[0]+"_"+selected_chain+'.pdb' # name output file
        with open (output_file, 'w') as pdb_splitted: # create for writing
            with open (current_file) as pdb_not_splitted: # open for reading
                for line in pdb_not_splitted:
                    if line.startswith('ATOM') and line[21]==selected_chain:
                        # select residues of the chain
                        pdb_splitted.write(line) # write in the output file
        sys.stderr.write("%s PDB file splitted successfully, chain %s selected \n" %(file_name[0], file_name[1]))
pdb_split_by_chain=select_chain_from_pdb(list_chain_separated, dir_pdb, dir_split)
