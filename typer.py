import pandas

def perType(filename):
    data = pandas.read_csv(filename, sep='\t', header=1)
    header = list(data.columns)

    # remove 'Position' and 'Reference' from column name list
    samples = header[2:]

    for sample in samples:
        print(list(data[sample]))




perType("variantContentTable.tsv")