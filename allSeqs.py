import sqlite3
import pandas
import os
import sqlparse


def createDB(db_name, tsvFile):
    # establish connection to sqlite database
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()

    if os.path.isfile(tsvFile):
        with open(tsvFile):
            # insert snp file into database
            pandas.read_csv(tsvFile, sep='\t').to_sql("mult_snps", conn, index=False, if_exists='replace')

    cursor.execute("ALTER TABLE mult_snps RENAME COLUMN 'SNP pattern' TO SNP_pattern;")

    sql_command = "SELECT * FROM mult_snps"
    # cursor.execute(sql_command)
    content = cursor.fetchall()
    print(content)

    conn.commit()

    conn.close()


# createDB("MultSnps.db", "allSeqsSnps")


def uniqSnpCount(db_name):
    # establish connection to sqlite database
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()

    sql_command = "SELECT * FROM mult_snps LIKE 'sequence%Contig';"

    cursor.execute(sql_command)
    content = cursor.fetchall()
    print(content)
    print(len(content))

    conn.commit()

    conn.close()

uniqSnpCount("Multsnps.db")
