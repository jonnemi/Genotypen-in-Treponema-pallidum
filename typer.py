import pandas
import os
import sqlite3

def uniqueSNPType(tsvFile):
    # establish connection to sqlite database
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()


    if os.path.isfile(tsvFile):
        with open(tsvFile):
            # insert snp file into database
            pandas.read_csv(tsvFile, sep='\t').to_sql("mult_snps", conn, index=False, if_exists='replace')


    # remove 'Position' and 'Reference' from column name list
    samples = header[2:]

    for sample in samples:
        print(list(data[sample]))




uniqueSNPType("A2_1.variants.tsv")