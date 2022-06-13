import sqlite3
import pandas
import os


def createDB(db_name, directory):
    # establish connection to sqlite database
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()

    # cursor.execute("DROP TABLE snps")

    sql_command = "CREATE TABLE IF NOT EXISTS snps ('SNP pattern' text, sequence_1_Contig text, sequence_1_PosInContg " \
                  "int, sequence_1_GenWidePos1 int, sequence_2_Contig text, sequence_2_PosInContg int, " \
                  "sequence_2_GenWidePos2 int); "
    cursor.execute(sql_command)

    # iterate over files in that directory
    for filename in os.listdir(directory):
        file = os.path.join(directory, filename)
        # checking if it is a file
        if os.path.isfile(file):
            with open(file) as tsvFile:
                # insert snp file into database
                pandas.read_csv(tsvFile, sep='\t').to_sql("snps", conn, index=False, if_exists='append')
    cursor.execute("ALTER TABLE snps RENAME COLUMN 'SNP pattern' TO SNP_pattern;")

    conn.commit()

    conn.close()

# createDB("snps.db", "snps")

def filterUnambiguous(from_db_name):
    # establish connection to sqlite database
    conn = sqlite3.connect(from_db_name)
    cursor = conn.cursor()

    # cursor.execute("DROP TABLE unambiguous")

    # create new table unambiguous
    sql_command = "CREATE TABLE IF NOT EXISTS unambiguous (SNP_pattern text, sequence_1_Contig text, " \
                  "sequence_2_GenWidePos2 int, count int); "
    cursor.execute(sql_command)

    # copy all unambiguous entries from table snps into new table
    sql_command = "INSERT INTO unambiguous SELECT SNP_pattern, sequence_1_Contig, sequence_2_GenWidePos2, COUNT(*) c " \
                  "FROM snps GROUP BY SNP_pattern, sequence_2_GenWidePos2 HAVING c = 1;"
    cursor.execute(sql_command)

    # delete last column count in table unambiguous
    sql_command = "ALTER TABLE unambiguous DROP COLUMN count;"
    cursor.execute(sql_command)

    conn.commit()

    conn.close()

# filterUnambiguous("snps.db")

def snpCount(db_name):
    # establish connection to sqlite database
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()

    # cursor.execute("DROP TABLE totalCount")

    # create new table
    sql_command = "CREATE TABLE IF NOT EXISTS totalCount (sequence_1_Contig text, total_count int); "
    cursor.execute(sql_command)

    sql_command = ("INSERT INTO totalCount SELECT sequence_1_Contig, COUNT(*) c FROM snps GROUP by sequence_1_Contig "
                   "HAVING c > 0 ORDER BY sequence_1_Contig;")
    cursor.execute(sql_command)

    conn.commit()

    conn.close()

# snpCount("snps.db")

def uniqSnpCount(db_name):
    # establish connection to sqlite database
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()

    # cursor.execute("DROP TABLE uniqCount")

    # create new table
    sql_command = "CREATE TABLE IF NOT EXISTS uniqCount (sequence_1_Contig text, uniq_count int); "
    cursor.execute(sql_command)

    sql_command = ("INSERT INTO uniqCount SELECT sequence_1_Contig, COUNT(*) c FROM unambiguous GROUP by "
                   "sequence_1_Contig HAVING c > 0 "
                   "ORDER BY sequence_1_Contig;")
    cursor.execute(sql_command)

    conn.commit()

    conn.close()

# uniqSnpCount("snps.db")

def totalVsUniqSnpCount(db_name):
    # establish connection to sqlite database
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()

    # cursor.execute("DROP TABLE totalVsUniqCount")

    # create new table allVSuniq
    sql_command = "CREATE TABLE IF NOT EXISTS totalVsUniqCount (sequence_1_Contig text, " \
                  "total_snp_count int, uniq_snp_count int); "
    cursor.execute(sql_command)

    # copy all counts into new table
    sql_command = "INSERT INTO totalVsUniqCount (sequence_1_Contig, total_snp_count) SELECT * FROM totalCount;"
    cursor.execute(sql_command)

    sql_command = "UPDATE totalVsUniqCount SET uniq_snp_count = uniqCount.uniq_count FROM uniqCount " \
                  "WHERE totalVsUniqCount.sequence_1_Contig = uniqCount.sequence_1_Contig;"
    cursor.execute(sql_command)

    sql_command = "SELECT * FROM totalVsUniqCount"
    cursor.execute(sql_command)
    content = cursor.fetchall()
    print(content)
    print(len(content))

    conn.commit()

    conn.close()

# totalVsUniqSnpCount("snps.db")


def ambiguousSnpCount(db_name):
    # establish connection to sqlite database
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()

    cursor.execute("DROP TABLE ambiguousCount")

    # create new table allVSuniq
    sql_command = "CREATE TABLE IF NOT EXISTS ambiguousCount (sequence_1_Contig text, " \
                  "total_count int, one_duplicate int, less_five_duplicates int, more_five_duplicates int); "
    cursor.execute(sql_command)

    sql_command = "INSERT INTO ambiguousCount (sequence_1_Contig, total_count) SELECT * FROM totalCount ORDER BY " \
                  "sequence_1_Contig "
    cursor.execute(sql_command)
#############################################hier weiter

    sql_command = "UPDATE ambiguousCount " \
                  "SET one_duplicate = SELECT sequence_1_Contig FROM snps COUNT(*) c GROUP BY SNP_pattern, sequence_2_GenWidePos2 HAVING c = 2 "\
                  "WHERE ambiguousCount.sequence_1_Contig = snps.sequence_1_Contig;"
    cursor.execute(sql_command)


    sql_command = "SELECT * FROM ambiguousCount"
    cursor.execute(sql_command)
    content = cursor.fetchall()
    print(content)
    print(len(content))

    conn.commit()

    conn.close()


# ambiguousSnpCount("snps.db")