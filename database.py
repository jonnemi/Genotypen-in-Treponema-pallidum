import sqlite3
import pandas
import os
import pandas as pd

def loadDB():
    if not os.path.isfile("snps.db"):
        # creating database for Reference-genome SNPs...
        print("Creating database for Reference-genome SNPs...")
        createDB("snps.db", "paarweise_MAUVE_SNPs")
        filterUnambiguous("snps.db")

def getGaps(filename):
    directory = "pairwise_gaps"
    file = os.path.join(directory, filename)
    # checking if it is a file
    if os.path.isfile(file):
        with open(file) as tsvFile:
            # insert snp file into database
            columns = ['Genome', 'Contig', 'Position_in_Contig', 'GenomeWide_Position',	'Length', 'Unnamed1', 'Unnamed2', 'Unnamed3', 'Unnamed4', 'Unnamed5', 'Unnamed6']
            df = pandas.read_csv(tsvFile, sep='\t', header=1, names=columns)

            df = df[df['Contig'] == '[1,1139633]']
            df = df[['GenomeWide_Position', 'Length']]
            df['Range'] = df[['GenomeWide_Position', 'Length']].to_numpy().tolist()
            df['Range'] = df['Range'].apply(lambda x: [x[0] + y for y in range(0, x[1])])
            flat_pos = [x for xs in df['Range'].tolist() for x in xs]

            gaps = pd.DataFrame(flat_pos, columns=['GenomeWide_Position'])
            gaps['Gaps_ahead'] = 1
            gaps['Gaps_ahead'] = gaps['Gaps_ahead'].cumsum()
            print(gaps)
        tsvFile.close()
    return gaps


def createDB(db_name, directory):
    # establish connection to sqlite database
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()

    #cursor.execute("DROP TABLE snps")

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
                df = pandas.read_csv(tsvFile, sep='\t')

                """# determine gap positions to remove
                gap_df = getGaps(filename)

                # remove gaps aka insertions
                df = df[~df["sequence_2_GenWidePos2"].isin(gap_df['GenomeWide_Position'].tolist())]

                # substract gaps from SNP positions
                df['Gap_count'] = df["sequence_2_GenWidePos2"].apply(lambda x: gap_df[gap_df['GenomeWide_Position'] <= x].tail(1).values)
                df['Gap_count'] = df['Gap_count'].apply(lambda x: x.tolist()[0][1] if x.tolist() else 0)
                df["sequence_2_GenWidePos2"] = df["sequence_2_GenWidePos2"] - df["Gap_count"]


                df.to_csv("test.tsv", sep='\t')
                df = df.drop(columns=["Gap_count"])"""

                df.to_sql("snps", conn, index=False, if_exists='append')
               

    # rename for easier access
    cursor.execute("ALTER TABLE snps RENAME COLUMN 'SNP pattern' TO SNP_pattern;")
    # rename columns for readability
    cursor.execute("ALTER TABLE snps RENAME COLUMN sequence_1_Contig TO sequence_name;")
    cursor.execute("ALTER TABLE snps RENAME COLUMN sequence_1_PosInContg TO sequence_PosInContg;")
    cursor.execute("ALTER TABLE snps RENAME COLUMN sequence_1_GenWidePos1 TO sequence_GenWidePos;")
    cursor.execute("ALTER TABLE snps RENAME COLUMN sequence_2_Contig TO reference_name;")
    cursor.execute("ALTER TABLE snps RENAME COLUMN sequence_2_PosInContg TO reference_PosInContg;")
    cursor.execute("ALTER TABLE snps RENAME COLUMN sequence_2_GenWidePos2 TO reference_GenWidePos;")

    conn.commit()

    conn.close()

#createDB("snps.db", "paarweise_MAUVE_SNPs")

def filterUnambiguous(from_db_name):
    # establish connection to sqlite database
    conn = sqlite3.connect(from_db_name)
    cursor = conn.cursor()

    #cursor.execute("DROP TABLE unambiguous")

    # create new table unambiguous
    sql_command = "CREATE TABLE IF NOT EXISTS unambiguous (SNP_pattern text, sequence_name text, " \
                  "reference_GenWidePos int); "
    cursor.execute(sql_command)

    # copy all unambiguous entries from table snps into new table
    sql_command = "INSERT INTO unambiguous SELECT SNP_pattern, sequence_name, reference_GenWidePos FROM snps;"
    cursor.execute(sql_command)

    ambigBases = ['r', 'y', 's', 'w', 'k', 'm', 'b', 'd', 'h', 'v', 'n']
    for base in ambigBases:
        sql_command = "DELETE FROM unambiguous WHERE SNP_pattern LIKE '%" + base + "%';"
        cursor.execute(sql_command)

    conn.commit()

    conn.close()

#filterUnambiguous("snps.db")

def filterUnique(from_db_name):
    # establish connection to sqlite database
    conn = sqlite3.connect(from_db_name)
    cursor = conn.cursor()

    cursor.execute("DROP TABLE unique")

    # create new table unique
    sql_command = "CREATE TABLE IF NOT EXISTS unique (SNP_pattern text, sequence_name text, " \
                  "reference_GenWidePos int, count int); "
    cursor.execute(sql_command)

    # copy all unique entries from table snps into new table
    sql_command = "INSERT INTO unique SELECT SNP_pattern, sequence_name, reference_GenWidePos, COUNT(*) c " \
                  "FROM snps GROUP BY SNP_pattern, reference_GenWidePos HAVING c = 1;"
    cursor.execute(sql_command)

    ambigBases = ['r', 'y', 's', 'w', 'k', 'm', 'b', 'd', 'h', 'v', 'n']
    for base in ambigBases:
        sql_command = "DELETE FROM unique WHERE SNP_pattern LIKE '%" + base + "%';"
        cursor.execute(sql_command)

    # delete last column count in table unique
    sql_command = "ALTER TABLE unique DROP COLUMN count;"
    cursor.execute(sql_command)

    conn.commit()

    conn.close()

# filterUnique("snps.db")

def snpCount(db_name):
    # establish connection to sqlite database
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()

    cursor.execute("DROP TABLE totalCount")

    # create new table
    sql_command = "CREATE TABLE IF NOT EXISTS totalCount (sequence_name text, total_count int); "
    cursor.execute(sql_command)

    sql_command = ("INSERT INTO totalCount SELECT sequence_name, COUNT(*) c FROM snps GROUP by sequence_name "
                   "HAVING c > 0 ORDER BY sequence_name;")
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
    sql_command = "CREATE TABLE IF NOT EXISTS uniqCount (sequence_name text, uniq_count int); "
    cursor.execute(sql_command)

    sql_command = ("INSERT INTO uniqCount SELECT sequence_name, COUNT(*) c FROM unique GROUP by "
                   "sequence_name HAVING c > 0 "
                   "ORDER BY sequence_name;")
    cursor.execute(sql_command)

    conn.commit()

    conn.close()

# uniqSnpCount("snps.db")

def totalVsUniqSnpCount(db_name):
    # establish connection to sqlite database
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()

    cursor.execute("DROP TABLE totalVsUniqCount")

    # create new table allVSuniq
    sql_command = "CREATE TABLE IF NOT EXISTS totalVsUniqCount (sequence_name text, " \
                  "total_snp_count int, uniq_snp_count int); "
    cursor.execute(sql_command)

    # copy all counts into new table
    sql_command = "INSERT INTO totalVsUniqCount (sequence_name, total_snp_count) SELECT * FROM totalCount;"
    cursor.execute(sql_command)

    sql_command = "UPDATE totalVsUniqCount SET uniq_snp_count = uniqCount.uniq_count FROM uniqCount " \
                  "WHERE totalVsUniqCount.sequence_name = uniqCount.sequence_name;"
    cursor.execute(sql_command)

    sql_command = "SELECT * FROM totalVsUniqCount"
    cursor.execute(sql_command)
    content = cursor.fetchall()
    print(content)
    print("Length of DB (should be 71): " + str(len(content)))

    sql_command = "SELECT * FROM totalVsUniqCount WHERE uniq_snp_count > 1"
    cursor.execute(sql_command)
    content = cursor.fetchall()
    print("Number of sequences with unique SNPs: " + str(len(content)))

    conn.commit()

    conn.close()

# totalVsUniqSnpCount("snps.db")