##########################################
############ PDB IDs in PFSAM ############
##########################################
### This program downloads family alignemnts from PFAM using its ID 

from urllib import request
import sys

def down_pfam(id_pfam):
    """ Downloads the alignment of a family of proteins from PFAM, given an id"""

    pfam_url = 'http://pfam.xfam.org/family/' + id_pfam + '/alignment/full'
    local_alig = 'alignment.sto'
    
    sys.stderr.write("Downloading the file...\n")
    request.urlretrieve(pfam_url,local_alig)
    
    # if the ID is not correct, an error file is going to be downloaded
    # check whether the file is the correct one or the error file
    
    with open (local_alig) as file_down:
        for line in file_down:
            if line.startswith("# STOCKHOLM 1.0"):
                sys.stderr.write("File download successfully \n")
                break # just look for the first line
                
            else:
                raise SystemExit("The PFAM ID was not found.")
                
