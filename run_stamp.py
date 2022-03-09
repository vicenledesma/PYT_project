######################################
############## RUN STAMP #############
######################################

import os

def run_stamp_roughfit (stamp_in_file, outFiles_prefix):

    """ Runs STAMP, a program for making structural multiple sequence 
    alignments, from the command line. It must be provided with an input
    file that has an understandable format for STAMP. Uses STAMP's ROUGHFIT
    option, which doesn't use any initial sequence information. The names of the output
    files are defined with the outFiles_prefix argument """

    # write command

    stamp_cmd = "stamp -l " + stamp_in_file + " -rough -n 2 -prefix " + outFiles_prefix + " > " + outFiles_prefix + ".out"    

    # execute command

    os.system(stamp_cmd)        

### EXAMPLE

run_stamp_roughfit("example.in", "stampoutfile")