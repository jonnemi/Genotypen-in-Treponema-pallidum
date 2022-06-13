import pandas


def readVCT(filename):
    data = pandas.read_csv(filename, sep='\t', header=1, index_col=False)
    print(data)


readVCT("variantContentTable.tsv")