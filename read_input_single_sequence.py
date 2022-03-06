######################################
############ READ SEQUENCE ###########
######################################

import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from urllib import request

user_input = sys.argv[1]
print(isinstance(user_input, str))

def read_input_seq(FASTA_or_ID):

    input_records = [] # empty array for storing input records

    ### Check the type of input
    ## The input is a FASTA file


    if isinstance(user_input, str):

        sys.stderr.write("UniProt ID provided, retrieving sequence as FASTA file...\n")

        uniprot_url = 'http://www.uniprot.org/uniprot/' + user_input + '.fasta'
        local_fasta = 'uniprot_seq.fasta'

        try:
            request.urlretrieve(uniprot_url, local_fasta)
        except:
            raise SystemExit('The ID provided could not be found in UniProt. Exiting the program.')

        user_input = 'uniprot_seq.fasta' # use the retrieved UniProt file as a FASTA input

    ## The input is a FASTA file

    elif os.path.isfile(user_input):

        sys.stderr.write("Reading protein sequence from FASTA file...\n")

        with open(user_input) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                input_records.append(record)



    ## Check the input for errors

   #if not input_records:
        #raise SystemExit('No input sequences found. Check that the sequences are in FASTA format. Exiting the program.')
    #elif len(input_records) > 1:
        #raise SystemExit('More than one input provided. Exiting the program.')
    #else:
        #return(input_records[0])

example = read_input_seq(user_input)
