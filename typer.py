import pandas
import os
import sqlite3


def uniqueSNPType(tsvFile):
    # establish connection to sqlite database
    conn = sqlite3.connect("snps.db")
    cursor = conn.cursor()

    # transform input file into pandas dataframe
    if os.path.isfile(tsvFile):
        query = pandas.read_csv(tsvFile, sep='\t', header=None, skiprows=1)
        query.columns = ['Position', 'SNP']

    # get all possible sequence_names aka types
    sql_command = "SELECT DISTINCT sequence_name FROM unambiguous;"
    cursor.execute(sql_command)
    content = cursor.fetchall()
    # list holding types and how often the query matched a SNP each type
    type_count = []
    for type in content:
        dic = {'type': type[0], 'count': 0}
        type_count.append(dic)

    query.reset_index()
    # print(listquery["SNP"].astype(str))
### das geht so nicht weil kombis durcheinander kommen
    sql_command = "SELECT sequence_name FROM unambiguous WHERE " \
                  "SNP_pattern IN " + str(tuple(query["SNP"])) + ";"
    print(sql_command)
    cursor.execute(sql_command)
    content = cursor.fetchall()
    print(content)





#uniqueSNPType("PT_SIF0908.variants.tsv")
