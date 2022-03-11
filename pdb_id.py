##########################################
############ PDB IDs in PFSAM ############
##########################################
### This program searchs por PDB IDs in the file downloaded from PFAM

from urllib import request
from Bio import SeqIO
from Bio import PDB
import sys

def get_id_pdb(alignment_file):
    
    """ Extract the PDB IDs and chains from the PFAM file. It stores this data 
    in a list of tupples, where the first position of the tuple is the PDB file 
    and the second one is the chain """
    
    list_id=[] # initialize list of tupples
    
    with open (alignment_file) as file_sto:
        for line in file_sto:
            if 'PDB' in line:
                id_tuple = (line[41:45], line[46:47])
                list_id.append(id_tuple) # append tupples to list
                                         # with ID and chain from PDB
    return (list_id)


def download_pdb(list_id_chain):
    
    """ Creates a set from the previous list, so the PDB IDs are unique
    (decrease computational cost). Then it downloads the PDB files for each of this IDs """
    
    set_unique_id = set() # avoid repeated PDB IDs

    for protein_id in list_id_chain:
        pdb_id = protein_id[0]
        set_unique_id.add(pdb_id)
        
    for id_unique in set_unique_id:
        PDB_file = '/Users/sarapoloalonso/Documents/Apuntes_2/PYT/Proyecto/pfam/pdb/' + id_unique + '.pdb'
        
        try:
            # Download pdb files
            request.urlretrieve('https://files.rcsb.org/download/' + id_unique + '.pdb', PDB_file)
            sys.stderr.write("%s PDB file downloaded successfully\n" %id_unique)
        
        except:
            raise SystemExit("ID PDB %s was not found." %id_unique)

            
def domains_templates(list_id_chain):
    
    """ Creates a text file with the needed format to be the input file
    of STAMP (structural alignment) """
    
    with open ("templates.domains", 'w') as stamp_input:
        for id_chain in list_id_chain:
            stamp_input.write('/pdb/'+id_chain[0]+'.pdb'+' '+id_chain[0]+' {CHAIN '+id_chain[1].upper()+'}\n')
            
