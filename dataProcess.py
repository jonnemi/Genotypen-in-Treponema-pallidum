import os
import sqlite3
import numpy as np
import pandas as pd
from Bio import SeqIO
import ast
from collections import Counter
pd.options.mode.chained_assignment = None  # default='warn'


def encodeNuc(nuc, encoding):
    symbols = ["a", "c", "g", "t", ".", "-"]
    i = symbols.index(nuc)
    if i == 5:
        i = 4
    if encoding == "string/none":
        enc = nuc
    elif encoding == "binary":
        if nuc != ".":
            enc = 1
        else:
            enc = 0
    elif encoding == "one-hot":
        one_hot_enc = [0] * 5
        one_hot_enc[i] = 1
        enc = one_hot_enc
    else:
        enc = i
    return enc

def getSNPdf(db_name):
    # get tSNE input data from snp database
    # establish connection to sqlite database
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()

    # get all possible sequence_names aka types
    sql_command = "SELECT DISTINCT sequence_name FROM unambiguous ORDER BY sequence_name;"
    cursor.execute(sql_command)
    content = cursor.fetchall()
    sequence_names = list()
    for seq in content:
        sequence_names.append(seq[0])

    # get all reference_GenWidePos for SNPs
    sql_command = "SELECT DISTINCT reference_GenWidePos FROM unambiguous ORDER BY reference_GenWidePos;"
    cursor.execute(sql_command)
    content = cursor.fetchall()
    positions = list()
    for pos in content:
        positions.append(pos[0])

    query = []
    for sequence in sequence_names:
        # query SNP pattern and position from database
        sql_command = "SELECT sequence_name, SNP_pattern, reference_GenWidePos FROM unambiguous WHERE" \
                      " unambiguous.sequence_name='" + str(sequence) + "' ORDER BY reference_GenWidePos;"
        cursor.execute(sql_command)
        content = cursor.fetchall()
        query.append(content)
    flat_query = [x for xs in query for x in xs]
    df = pd.DataFrame(flat_query, columns=['Sequence_name', 'SNP', 'Position'])
    df['Reference'] = df['SNP'].str[1]
    df['SNP'] = df['SNP'].str[0]
    return df


def encDF(df, encoding, common_pos):
    if encoding == "binary":
        default_enc = 0
        enc_len = 1
    elif encoding == "one-hot":
        default_enc = [0, 0, 0, 0, 1]
        enc_len = len(default_enc)
    elif encoding == "integer":
        default_enc = 4
        enc_len = 1
    elif encoding == "string/none":
        default_enc = "."
        enc_len = len(default_enc)
    else:
        raise ValueError('Please specify encoding: (binary, one-hot, integer, string/none')

    seq_count = df['Sequence_name'].nunique()
    sequences = df['Sequence_name'].unique().tolist()

    # create training vector of SNPs, with field for each SNP site (reference_GenWidePos) using one-hot encoding
    reference_3D = np.full((1, len(common_pos), enc_len), default_enc)
    vec_3D = np.full((seq_count, len(common_pos), enc_len), default_enc)

    for index, row in df.iterrows():
        if row['Position'] in common_pos:
            seq = sequences.index(row["Sequence_name"])
            pos = common_pos.index(row["Position"])

            enc_entry = encodeNuc(row["SNP"], encoding)
            vec_3D[seq][pos] = enc_entry

            reference_3D[0][pos] = default_enc #encodeNuc(row["Reference"], encoding)

    vec_3D = np.append(vec_3D, reference_3D, axis=0)

    # reshape training vector from 3D to 2D
    sample, position, one_hot_enc = vec_3D.shape
    vec_2D = vec_3D.reshape((sample, position * one_hot_enc))
    return vec_2D


def encQuery(df, encoding, common_pos):
    if encoding == "binary":
        default_enc = 0
        enc_len = 1
    elif encoding == "one-hot":
        default_enc = [0, 0, 0, 0, 1]
        enc_len = len(default_enc)
    elif encoding == "integer":
        default_enc = 4
        enc_len = 1
    elif encoding == "string/none":
        default_enc = "."
        enc_len = len(default_enc)
    else:
        raise ValueError('Please specify encoding: (binary, one-hot, integer, string/none')

    # drop all columns aka samples with no SNPs
    df = df.drop(columns=df.columns[(df == '.').all()])

    sample_names = list(df.columns)
    sample_names.remove('Position')
    sample_names.remove('Reference')

    # create vector query SNPs, with field for each SNP site (reference_GenWidePos) using one-hot encoding
    eval_3D = np.full((len(sample_names) + 1, len(common_pos), enc_len), default_enc)

    s = 0
    for sample in sample_names:
        p = 0
        for entry in df[sample]:
            pos = common_pos.index(df['Position'].iloc[p])
            eval_3D[s][pos] = encodeNuc(entry.lower(), encoding)
            p += 1
        s += 1

    # add vector for Reference containing only default symbols
    for c in range(len(common_pos) - 1):
        eval_3D[s][c] = default_enc
    sample_names.append('Reference')

    # reshape training vector from 3D to 2D
    sample, position, one_hot_enc = eval_3D.shape
    eval_2D = eval_3D.reshape((sample, position * one_hot_enc))
    query = {'data': eval_2D, 'sample_names': sample_names}
    return query


def newEncQuery(df, encoding, common_pos):
    if encoding == "binary":
        default_enc = 0
        enc_len = 1
    elif encoding == "one-hot":
        default_enc = [0, 0, 0, 0, 1]
        enc_len = len(default_enc)
        common_pos = list(set(common_pos))
    elif encoding == "integer":
        default_enc = 4
        enc_len = 1
        common_pos = list(set(common_pos))
    elif encoding == "string/none":
        default_enc = "."
        enc_len = len(default_enc)
        common_pos = list(set(common_pos))
    else:
        raise ValueError('Please specify encoding: (binary, one-hot, integer, string/none')

    sample_names = df["SAMPLES"].str.split(",").tolist()
    sample_names = [item for sublist in sample_names for item in sublist]
    sample_names = list(set(sample_names))
    sample_names.sort()

    position_len = len(common_pos)

    # create vector query SNPs, with field for each SNP site (reference_GenWidePos)
    eval_3D = np.full((len(sample_names) + 1, position_len, enc_len), default_enc)
    f = []
    for position in df['Position']:
        # only keep SNPs of impact given in filter (LOW, MODERATE or HIGH)
        pos = common_pos.index(position)
        if encoding == "binary":
            common_pos[pos] = -1
        sequences = df[df['Position'] == position]["SAMPLES"].str.split(",").tolist()
        sequences = sequences[0]

        for sequence in sequences:
            seq = sample_names.index(sequence)
            if encoding == "binary":
                eval_3D[seq][pos] = 1
            else:
                entry = df[df['Position'] == position]['ALT_CONTENT'].values[0]
                if '-' not in entry and len(entry) == 1:
                    # only consider SNPs, no indels
                    eval_3D[seq][pos] = encodeNuc(entry.lower(), encoding)
        pos += 1

    # add vector for Reference containing only default symbols
    for c in range(position_len - 1):
        eval_3D[len(sample_names)][c] = default_enc
    sample_names.append('Reference')

    # reshape training vector from 3D to 2D
    sample, position, one_hot_enc = eval_3D.shape
    eval_2D = eval_3D.reshape((sample, position * one_hot_enc))
    query = {'data': eval_2D, 'sample_names': sample_names}
    return query


def assignStrains(sample_names):
    # all samples manually grouped by strain
    group1 = ['BosniaA']
    group2 = ['Fribourg', 'CDC2', 'GHA1', 'Gauthier', 'IND1', 'SAM1', 'SamoaD']
    group3 = ['Seattle81', 'SEA86', 'NE20', 'BAL3', 'BAL73', 'Chicago', 'NIC1', 'NIC2', 'Dallas']

    strain_labels = ['TEN', 'TPE', 'Nichols-like', 'SS14-like']

    strains = []
    for sample in sample_names:
        if sample in group1:
            strain = strain_labels[0]
        elif sample in group2:
            strain = strain_labels[1]
        elif sample in group3:
            strain = strain_labels[2]
        else:
            strain = strain_labels[3]
        strains.append(strain)
    return strains


def getADAdataset(tsvFile, db_name, encoding):
    # get all SNPs from reference database
    SNPdf = getSNPdf(db_name)
    ref_positions = SNPdf["Position"].unique().tolist()
    ref_positions.sort()

    # read tsvFile containing query sequence SNPs
    if os.path.isfile(tsvFile):
        query = pd.read_csv(tsvFile, sep='\t', header=1)

    # manually drop TPE an TEN samples
    query = query.drop(columns=['BosniaA', 'Fribourg', 'CDC2', 'GHA1', 'Gauthier', 'IND1', 'SAM1', 'SamoaD'])

    # create vector of SNPs, with field for each SNP site (reference_GenWidePos)
    sample_names = list(query.columns)
    sample_names.remove('Position')
    sample_names.remove('Reference')

    # find common SNP positions in reference database and query
    query_pos = query.Position.tolist()
    common_pos = list(set(map(str, ref_positions)).intersection(query_pos))
    common_pos.sort(key=lambda pos: int(pos))

    # remove all rows from dataframes that don't have common reference genome position for SNP
    query = query[query['Position'].isin(common_pos)]
    query['Position'] = pd.to_numeric(query['Position'])
    query.reset_index()
    common_pos = list(map(int, common_pos))
    SNPdf = SNPdf[SNPdf['Position'].isin(common_pos)]
    SNPdf.reset_index()

    sequence_names = SNPdf['Sequence_name'].unique().tolist()
    sequence_names.sort()
    sequence_names.append("Reference(NC_021490)")

    # create training and test vector of SNPs, with field for each SNP site (reference_GenWidePos) using one-hot encoding
    train_2D = encDF(SNPdf, encoding, common_pos)

    test_process = encQuery(query, encoding, common_pos)
    test_2D = test_process["data"]
    names = test_process["sample_names"]

    dataset = {'train': train_2D, 'train_seq_names': sequence_names, 'test': test_2D, 'test_seq_names': names}
    return dataset


def getLociDataset(tsvFile, db_name, default_enc, loci_list, new_format, filter=[]):
    # load loci dataframe
    loci_df = pd.read_csv("loci.tsv", sep='\t', converters={"var sites": ast.literal_eval})

    # create list of MLST loci SNP positions
    MLST_loci = loci_df[loci_df['locus_tag'].isin(loci_list)]
    MLST_positions = []
    MLST_loci["var sites"].apply(lambda x: MLST_positions.append(list(range(min(x) + 1, max(x) + 1))))

    # manually add 23S rRNA gene positions (gene is duplicated), because locus isn't listed in loci.tsv
    rRNA_start_1 = 233099
    rRNA_end_1 = 236051
    rRNA_start_2 = 281544
    rRNA_end_2 = 284496


    rRNA_pos = list(range(rRNA_start_1, rRNA_end_1 + 1))
    rRNA_pos.extend(list(range(rRNA_start_2, rRNA_end_2 + 1)))

    MLST_positions.append(rRNA_pos)
    for x in MLST_positions:
        x.sort()

    flat_MLST = [item for sublist in MLST_positions for item in sublist]

    # get all SNPs from reference database
    SNPdf = getSNPdf(db_name)
    #DBdf_toMega(SNPdf, flat_MLST)

    # read tsvFile containing query sequence SNPs
    if not new_format:
        if os.path.isfile(tsvFile):
            query = pd.read_csv(tsvFile, sep='\t', header=1)

        # manually drop TPE an TEN samples
        query = query.drop(columns=['BosniaA', 'Fribourg', 'CDC2', 'GHA1', 'Gauthier', 'IND1', 'SAM1', 'SamoaD'])
        position_col = 'Position'

        query = query.replace("-", ".")

        all_sample_names = list(query.columns)
        all_sample_names.remove(position_col)


        # print query SNPs to .meg file
        #oldQuerydf_toMega(query, flat_MLST)
    else:
        if os.path.isfile(tsvFile):
            query = pd.read_csv(tsvFile, sep='\t', header=0)
        query.rename(columns = {'POSITION':'Position'}, inplace = True)

        # print query SNPs to .meg file
        #newQuerydf_toMega(query, flat_MLST)

        # remove all indel rows from dataframe
        query = query[query['ALT_CONTENT'].str.len() == 1 & ~query['ALT_CONTENT'].str.contains("-")]
        query.reset_index()

        # filter for impact
        if filter:
            query['IMPACT'] = query["INFO"].apply(lambda x: x.split("IMPACT=")[1].split(",")[0])
            query = query[query["IMPACT"].isin(filter)].reset_index()

        query['Position'] = query['Position'].astype(str)
        all_sample_names = query["SAMPLES"].apply(lambda x: x.split(",")).tolist()
        all_sample_names = [item for sublist in all_sample_names for item in sublist]
        all_sample_names = list(set(all_sample_names))
        all_sample_names.append("Reference")


    sample_names = []

    train_loci = []
    test_loci = []
    train_names = []
    # create vector of SNPs, with field for each locus SNP site
    for locus in MLST_positions:
        # remove all rows from dataframes that don't have common reference genome position for SNP
        # samples where the whole locus is missing will later have an X in the allelic profile
        locus_query = query[query['Position'].isin(map(str, locus))]
        locus_query['Position'] = pd.to_numeric(locus_query['Position'])
        locus_query.reset_index()

        locus_SNPdf = SNPdf[SNPdf['Position'].isin(locus)]
        locus_SNPdf.reset_index()

        sequence_names = locus_SNPdf['Sequence_name'].unique().tolist()
        sequence_names.sort()
        sequence_names.append("Reference")
        train_names.append(sequence_names)

        # create training and test vector of SNPs, with field for each SNP site (reference_GenWidePos) using one-hot encoding
        train_2D = encDF(locus_SNPdf, default_enc, locus)
        if new_format:
            test_process = newEncQuery(locus_query, default_enc, locus)
        else:
            test_process = encQuery(locus_query, default_enc, locus)

        test_2D = test_process["data"]
        names = test_process["sample_names"]
        sample_names.append(names)

        train_loci.append(train_2D)
        test_loci.append(test_2D)

    loci_list.append("23S rRNA")
    all_train_names = SNPdf['Sequence_name'].unique().tolist()
    all_train_names.append("Reference")
    dataset = {'train': train_loci, 'train_seq_names': train_names, 'all_train_names': all_train_names,
               'test': test_loci, 'test_seq_names': sample_names, 'all_test_names': all_sample_names, 'loci': loci_list}
    return dataset

def getStrains(id_list):
    directory = "RefSeqGBKs"
    strains = []
    for filename in os.listdir(directory):
        if filename.split('.')[0] in id_list or "Reference" in filename:
            file = os.path.join(directory, filename)
            # checking if it is a file
            if os.path.isfile(file):
                try:
                    genome_record = SeqIO.read(file, "genbank")
                    for feature in genome_record.features:
                        if feature.type == "source":
                            strain = feature.qualifiers.get("strain")
                            strains.append(strain[0])
                except:
                    print("Reading error in " + file + ", please look up strain manually.")
    return strains

# print(encDF(getSNPdf("snps.db"), [0, 0, 0, 0, 1], com_pos))
# getLociDataset("66-Proben-Datensatz.tsv", "snps.db", [0, 0, 0, 0, 1])


def getQueryDataset(tsvFile, filter, encoding):
    if os.path.isfile(tsvFile):
        query = pd.read_csv(tsvFile, sep='\t', header=0)
    query.rename(columns={'POSITION': 'Position'}, inplace=True)

    # create vector of SNPs, with field for each SNP site
    # remove all indel rows from dataframe
    query = query[query['ALT_CONTENT'].str.len() == 1 & ~query['ALT_CONTENT'].str.contains("-")]
    query.reset_index()

    # filter for impact
    query['IMPACT'] = query["INFO"].apply(lambda x: x.split("IMPACT=")[1].split(",")[0])
    query = query[query["IMPACT"].isin(filter)].reset_index()

    """df = pd.DataFrame(query)
    df.to_csv("query.tsv", sep='\t')"""

    return query

def checkContain(samples, mlst_samples):
    contains = False
    for sample in samples:
        if sample in mlst_samples:
            contains = True
    return contains


def newQuerydf_toMega(df, MLST_positions, file="newQuery_full.meg"):
    # filter out all irrelevant samples for MLST
    tsvFile = "query_MLST.csv"
    if os.path.isfile(tsvFile):
        mlst = pd.read_csv(tsvFile, sep=',', header=0)
    mlst_samples = ["Reference"]
    mlst['Samples'].apply(lambda x: mlst_samples.extend(x.split(", ")))

    df = df[df['ALT_CONTENT'].str.len() == 1 & ~df['ALT_CONTENT'].str.contains("-")]
    df = df[df['Position'].isin(MLST_positions)].reset_index(drop=True)
    df['SAMPLES'] = df['SAMPLES'].apply(lambda x: x.split(","))
    df['MLST_SAMPLES'] = df['SAMPLES'].apply(lambda x: checkContain(x, mlst_samples))
    df = df[df['MLST_SAMPLES'] == True]

    positions = df["Position"].unique().tolist()

    strings = np.full((len(mlst_samples), len(positions)), ".")

    # print SNP sequences to MEGA file
    with open(file, 'w') as f:
        print("#mega", file=f)
        print("TITLE: Noninterleaved sequence data", file=f)
        print(file=f)

        # iterate over all positions, hence all SNPs
        for index, row in df.iterrows():
            pos = row['Position']
            pos_idx = positions.index(pos)
            alt_base = row['ALT_CONTENT']

            ref = row['REF_CONTENT']
            seqs = row['SAMPLES']

            for seq in seqs:
                if seq in mlst_samples:
                    seq_idx = mlst_samples.index(seq)
                    strings[0][pos_idx] = ref
                    strings[seq_idx][pos_idx] = alt_base

        for sample in mlst_samples:
            idx = mlst_samples.index(sample)
            label = "#" + sample
            print(label, file=f)

            sampleSNPs = "".join(strings[idx])
            #sampleSNPs = sampleSNPs.replace("-", ".")
            print(sampleSNPs, file=f)

    print("Mega file of SNP sequences was saved to " + file)


def oldQuerydf_toMega(df, MLST_positions, file="oldQuery_full.meg"):

    df = df[df['Position'].isin(map(str, MLST_positions))].reset_index(drop=True)

    samples = list(df.columns)
    samples.remove('Position')

    # print SNP sequences to MEGA file
    with open(file, 'w') as f:
        print("#mega", file=f)
        print("TITLE: Noninterleaved sequence data", file=f)
        print(file=f)

        # iterate over all positions, hence all SNPs
        for sample in samples:
            label = "#" + sample
            print(label, file=f)

            sampleSNPs = "".join(df[sample].tolist())
            sampleSNPs = sampleSNPs.replace("-", ".")
            print(sampleSNPs, file=f)

    print("Mega file of SNP sequences was saved to " + file)


def DBdf_toMega(df, MLST_positions, file="refDB_full.meg"):
    samples = df['Sequence_name'].unique().tolist()
    samples.insert(0, "Reference")

    df = df[df['Position'].isin(MLST_positions)].reset_index(drop=True)
    positions = df['Position'].unique().tolist()
    positions.sort()

    # print SNP sequences to MEGA file
    with open(file, 'w') as f:
        print("#mega", file=f)
        print("TITLE: Noninterleaved sequence data", file=f)
        print(file=f)

        strings = np.full((len(samples), len(positions)), ".")

        # iterate over all positions, hence all SNPs
        for index, row in df.iterrows():
            pos = row['Position']
            pos_idx = positions.index(pos)
            alt_base = row['SNP']

            ref = row['Reference']
            seq_idx = samples.index(row['Sequence_name'])

            strings[0][pos_idx] = ref
            strings[seq_idx][pos_idx] = alt_base


        for sample in samples:
            idx = samples.index(sample)
            label = "#" + sample
            print(label, file=f)

            sampleSNPs = "".join(strings[idx])
            sampleSNPs = sampleSNPs.replace("-", ".")
            print(sampleSNPs, file=f)

    print("Mega file of SNP sequences was saved to " + file)

#getLociDataset("66-Proben-Datensatz.tsv", "snps.db", "string/none", ["TPANIC_RS00695", "TPANIC_RS02695", "TPANIC_RS03500"], False)
#getLociDataset("1508-Proben-Datensatz.tsv", "snps.db", "string/none", ["TPANIC_RS00695", "TPANIC_RS02695", "TPANIC_RS03500"], True)


def groupMetadata(tsvFile):
    # read metadata
    if os.path.isfile(tsvFile):
        data = pd.read_csv(tsvFile, sep='\t', header=0)
    # groupby country and institute
    group = data.groupby(['country', 'institute'])['strain'].size().reset_index(name='counts')

    # print to file
    print(group)
    dir = tsvFile.split("/")[0]
    filename = "grouped_" + tsvFile.split("/")[1]
    new_file = dir + "/" + filename
    print(new_file)
    group.to_csv(new_file, sep='\t', index=False)

#groupMetadata("metadata/62_TPA_Meta.tsv")
#groupMetadata("metadata/U19_TPA_Meta.tsv")


