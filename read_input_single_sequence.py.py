######################################
############ READ SEQUENCE ###########
######################################

import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from urllib import request

user_input = sys.argv[1]

def uniprot_retrieve(ID):

    """ Function that takes a UniProt ID and downloads its sequence in FASTA format
        from the online server"""

    sys.stderr.write("UniProt ID provided, retrieving sequence as FASTA file...\n")

    uniprot_url = 'http://www.uniprot.org/uniprot/' + user_input + '.fasta'
    local_fasta = 'uniprot_seq.fasta'

    try: # check if ID exists
        request.urlretrieve(uniprot_url, local_fasta)
    except:
        raise SystemExit('The ID provided could not be found in UniProt. Exiting the program.')

    return(local_fasta)

def FASTA_parse(FASTA_file):

    """ Function that parses a FASTA file and returns a Bio Record object """

    input_records = [] # empty array for storing input records

    sys.stderr.write("Reading protein sequence from FASTA file...\n")

    with open(FASTA_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            input_records.append(record)

    ## Check the input for errors

    if not input_records:
        raise SystemExit('No input sequences found. Check that the sequences are in FASTA format. Exiting the program.')
    elif len(input_records) > 1:
        raise SystemExit('More than one input provided. Exiting the program.')
    else:
        sys.stderr.write("Input processed corectly!\n")
        return(input_records[0])
