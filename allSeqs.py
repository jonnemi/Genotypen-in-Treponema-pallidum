import sqlite3
import pandas
import os


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

    # get all column names
    sql_command = "PRAGMA table_info(mult_snps);"

    cursor.execute(sql_command)
    content = cursor.fetchall()

    # filter genome wide positions for all sequences
    seqsGenWidePos = list()
    count = 0
    for columnInfo in content:
        if count % 3 == 0:
            seqsGenWidePos.append(columnInfo[1])
        count += 1

    # remove first entry, since it's not a genome wide position
    seqsGenWidePos.remove('SNP_pattern')

    sql_command = "SELECT " + ', '.join(seqsGenWidePos) + " FROM mult_snps;"
    cursor.execute(sql_command)
    content = cursor.fetchall()

    # count rows aka snps that are unique for all sequences concerning the genome wide position
    """ uniqueCount = 0
    for row in content:
        row = list(row)
        if 0 in row:
            row = [i for i in row if i != 0]
            row.append('dummy')
        if (len(set(row)) == len(row)):
            uniqueCount += 1"""

    uniqueCount = 0
    for row in content:
        row = set([x for x in row if row.count(x) > 1])
        row.discard(0)
        if len(row) == 0:
            uniqueCount += 1


    # get total snps count
    sql_command = "SELECT SNP_pattern FROM mult_snps;"
    cursor.execute(sql_command)
    content = cursor.fetchall()
    totalCount = len(content)

    print("Total number of SNPs from MSA: " + str(totalCount) + ", unique SNPs: " + str(uniqueCount))

    conn.commit()

    conn.close()

# uniqSnpCount("Multsnps.db")


def totalVsUniqSnpCount(db_name):
    # establish connection to sqlite database
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()

    # cursor.execute("DROP TABLE totalVsUniqCount")

    """# create new table allVSuniq
    sql_command = "CREATE TABLE IF NOT EXISTS totalVsUniqCount (sequence_1_Contig text, " \
                  "total_snp_count int, uniq_snp_count int); "
    cursor.execute(sql_command)

    # copy all counts into new table
    sql_command = "INSERT INTO totalVsUniqCount (sequence_1_Contig, total_snp_count) SELECT * FROM totalCount;"
    cursor.execute(sql_command)

    sql_command = "UPDATE totalVsUniqCount SET uniq_snp_count = uniqCount.uniq_count FROM uniqCount " \
                  "WHERE totalVsUniqCount.sequence_1_Contig = uniqCount.sequence_1_Contig;"
    cursor.execute(sql_command)"""

    sql_command = "SELECT * FROM mult_snps"
    cursor.execute(sql_command)
    content = cursor.fetchall()
    print(content)
    print(len(content))

    conn.commit()

    conn.close()

#totalVsUniqSnpCount("Multsnps.db")
