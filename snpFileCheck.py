import pandas
import os

def checkSNPFile(filename):
    my_data = pandas.read_csv(filename, sep='\t', index_col=False)
    if my_data["sequence_2_Contig"][0] != ("NC_021490"):
        raise Exception("Wrong reference genome position in " + filename)
    if not(my_data["sequence_1_GenWidePos1"].equals(my_data["sequence_1_PosInContg"])):
        raise Exception("Position deviation in seq in " + filename)
    if not(my_data["sequence_2_GenWidePos2"].equals(my_data["sequence_2_PosInContg"])):
        raise Exception("Position deviation in reference in " + filename)
    # print("SNP file " + filename  + "is fine.")
    return True


def checkSNPDir(directory):
    approved = False
    # iterate over files in that directory
    for filename in os.listdir(directory):
        file = os.path.join(directory, filename)
        # checking if it is a file
        if os.path.isfile(file):
            approved = checkSNPFile(file)
            if not approved:
                break
    return approved

# Usage:
# assign directory
print(checkSNPDir("snps"))

