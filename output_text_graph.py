###################################################
############## TEXT AND PLOT OUTPUT ###############
###################################################
# Libraries
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

def create_output_text (prefix, complete_flex_list):

    '''Creates parseable text file with the flexibility scores'''

    fd = open(prefix + ".txt", "w")
    fd.write("Position\t\t\tResidue\t\t\tFlexibility\t\tConfidence\n")

    for flex_tuple in complete_flex_list:
        fd.write(str(flex_tuple[0]) + "\t\t\t\t" + flex_tuple[1] + "\t\t\t" + "%+.3f\t\t\t" % (flex_tuple[2]) + str(flex_tuple[3]) + "\n")

    fd.close()


def draw_flex_line(prefix, flexibility_results):

    '''Creates plot of the flexibility results'''

    pos = []
    res = []
    flex = []
    conf = []

    for flex_tuple in flexibility_results:

        pos.append(flex_tuple[0])
        res.append(flex_tuple[1])
        flex.append(flex_tuple[2])
        conf.append(flex_tuple[3])

    plt.plot(pos,flex)
    plt.axhline(y=0, color='r', linestyle='-')
    plt.xlabel("Residue position")
    plt.ylabel("Flexibility score")

    plt.savefig(prefix + '.png')
    plt.show()

