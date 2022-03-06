from urllib import request
import sys
########## This program downloads alignemnt from PFAM of a pfam family ##########

####### The input is the ID of the pfam family ######

id_pfam=sys.argv[1]

def down_pfam(id_pfam):
    """Function that downloads the alignment of a family
    of proteins from PFAM, given an id"""

    pfam_url='http://pfam.xfam.org/family/'+id_pfam+'/alignment/full'
    local_alig = 'alignment.sto'
    sys.stderr.write("Downloading the file...\n")
    request.urlretrieve(pfam_url,local_alig)
    ##### if there PF is not correct, an error file is going to be downloaded
    ##### It is necessary to check whether the file is the correct one or not
    with open (local_alig) as file_down:
        for line in file_down:
            if line.startswith("# STOCKHOLM 1.0"):
                sys.stderr.write("File download successfully \n")
                break # just look for the first line
            else:
                raise SystemExit("The PFAM ID was not found.")
sara=down_pfam(id_pfam)
