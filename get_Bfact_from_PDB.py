################################################
############ GET B-FACTORS FROM PDB ############
################################################

from Bio.PDB import Atom
from Bio.PDB import PDBParser
import numpy

def get_B_from_PDB (PDB_file):

    """ Reads a PDB file and extracts the B factors for
    the alpha carbons of each amino acid. Returns the B factors
    as a list """

    # parse file

    parser = PDBParser()
    pdb_struct = parser.get_structure("wtf" ,PDB_file)

    # access atoms
    # create B factors list for alpha carbons

    B_list = []

    for model in pdb_struct:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.get_name() == "CA": # only for alpha carbons
                        B_list.append(atom.get_bfactor())
    
    return(B_list)