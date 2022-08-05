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
        return nuc
    elif encoding == "binary":
        return 1
    elif encoding == "one-hot":
        one_hot_enc = [0] * 5
        one_hot_enc[i] = 1
        return one_hot_enc
    else:
        return i

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


def encDF(df, default_enc, common_pos):
    if isinstance(default_enc, list):
        symb_count = len(default_enc)
    else:
        symb_count = 1

    seq_count = df['Sequence_name'].nunique()
    sequences = df['Sequence_name'].unique().tolist()

    # create training vector of SNPs, with field for each SNP site (reference_GenWidePos) using one-hot encoding
    reference_3D = np.full((1, len(common_pos), symb_count), default_enc)
    vec_3D = np.full((seq_count, len(common_pos), symb_count), default_enc)

    for index, row in df.iterrows():
        if row['Position'] in common_pos:
            seq = sequences.index(row["Sequence_name"])
            pos = common_pos.index(row["Position"])

            enc_entry = encodeNuc(row["SNP"], default_enc)
            vec_3D[seq][pos] = enc_entry

            reference_3D[0][pos] = encodeNuc(row["Reference"], default_enc)

    vec_3D = np.append(vec_3D, reference_3D, axis=0)

    # reshape training vector from 3D to 2D
    sample, position, one_hot_enc = vec_3D.shape
    vec_2D = vec_3D.reshape((sample, position * one_hot_enc))
    return vec_2D


def encQuery(df, default_enc, common_pos):
    if isinstance(default_enc, list):
        symb_count = len(default_enc)
    else:
        symb_count = 1

    sample_names = list(df.columns)
    sample_names.remove('Position')
    sample_names.remove('Reference')

    # create vector query SNPs, with field for each SNP site (reference_GenWidePos) using one-hot encoding
    eval_3D = np.full((len(sample_names), len(common_pos), symb_count), default_enc)

    s = 0
    for sample in sample_names:
        p = 0
        for entry in df[sample]:
            pos = common_pos.index(df['Position'].iloc[p])
            eval_3D[s][pos] = encodeNuc(entry.lower(), default_enc)
            p += 1
        s += 1

    # reshape training vector from 3D to 2D
    sample, position, one_hot_enc = eval_3D.shape
    eval_2D = eval_3D.reshape((sample, position * one_hot_enc))
    return eval_2D


def newEncQuery(df, encoding, common_pos):
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

    sample_names = df["SAMPLES"].str.split(",").tolist()
    sample_names = [item for sublist in sample_names for item in sublist]
    sample_names = list(set(sample_names))
    sample_names.sort()

    # create vector query SNPs, with field for each SNP site (reference_GenWidePos)
    eval_3D = np.full((len(sample_names), len(common_pos), enc_len), default_enc)

    pos = 0
    for position in df['Position']:
        # only keep SNPs of impact given in filter (LOW, MODERATE or HIGH)
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


def getADAdataset(tsvFile, db_name, default_enc):
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
    train_2D = encDF(SNPdf, default_enc, common_pos)
    test_2D = encQuery(query, default_enc, common_pos)
    dataset = {'train': train_2D, 'train_seq_names': sequence_names, 'test': test_2D, 'test_sample_names': sample_names}
    return dataset


def getLociDataset(tsvFile, db_name, default_enc, loci_list, new_format, filter):
    # load loci dataframe
    loci_df = pd.read_csv("loci.tsv", sep='\t', converters={"var sites": ast.literal_eval})
    # create list of MLST loci SNP positions
    MLST_loci = loci_df[loci_df['locus_tag'].isin(loci_list)]
    MLST_positions = []
    MLST_loci["var sites"].apply(lambda x: MLST_positions.append(x))

    rRNA_R8_v1 = 235204
    rRNA_R9_v1 = 235205
    rRNA_R8_v2 = 283649
    rRNA_R9_v2 = 283650

    rRNA_pos = [rRNA_R8_v1, rRNA_R9_v1, rRNA_R8_v2, rRNA_R9_v2]

    MLST_positions.append(rRNA_pos)
    for x in MLST_positions:
        x.sort()

    # get all SNPs from reference database
    SNPdf = getSNPdf(db_name)

    # read tsvFile containing query sequence SNPs
    if not new_format:
        if os.path.isfile(tsvFile):
            query = pd.read_csv(tsvFile, sep='\t', header=1)

        # manually drop TPE an TEN samples
        query = query.drop(columns=['BosniaA', 'Fribourg', 'CDC2', 'GHA1', 'Gauthier', 'IND1', 'SAM1', 'SamoaD'])
        position_col = 'Position'

        sample_names = list(query.columns)
        sample_names.remove(position_col)
        sample_names.remove('Reference')
    else:
        if os.path.isfile(tsvFile):
            query = pd.read_csv(tsvFile, sep='\t', header=0)
        query.rename(columns = {'POSITION':'Position'}, inplace = True)
        query['Position'] = query['Position'].astype(str)
        all_sample_names = query["SAMPLES"].apply(lambda x: x.split(",")).tolist()
        all_sample_names = [item for sublist in all_sample_names for item in sublist]
        all_sample_names = list(set(all_sample_names))
        all_sample_names.sort()

        sample_names = []

    train_loci = []
    test_loci = []
    train_names = []
    # create vector of SNPs, with field for each locus SNP site
    for locus in MLST_positions:
        # remove all rows from dataframes that don't have common reference genome position for SNP
        locus_query = query[query['Position'].isin(map(str, locus))]
        locus_query['Position'] = pd.to_numeric(locus_query['Position'])
        locus_query.reset_index()

        locus_SNPdf = SNPdf[SNPdf['Position'].isin(locus)]
        locus_SNPdf.reset_index()

        sequence_names = locus_SNPdf['Sequence_name'].unique().tolist()
        sequence_names.sort()
        sequence_names.append("Reference(NC_021490)")
        train_names.append(sequence_names)

        # create training and test vector of SNPs, with field for each SNP site (reference_GenWidePos) using one-hot encoding
        train_2D = encDF(locus_SNPdf, default_enc, locus)
        if new_format:
            test_process = newEncQuery(locus_query, default_enc, locus, filter)
            test_2D = test_process["data"]
            names = test_process["sample_names"]
            sample_names.append(names)
        else:
            test_2D = encQuery(locus_query, default_enc, locus)
        train_loci.append(train_2D)
        test_loci.append(test_2D)

    loci_list.append("23S rRNA")
    all_train_names = SNPdf['Sequence_name'].unique().tolist()
    all_train_names.append("Reference(NC_021490)")
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
# getLociDataset("variantContentTable.tsv", "snps.db", [0, 0, 0, 0, 1])


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


"""def SNPdf_toMega(df, path):
    sample_names = df["SAMPLES"].str.split(",").tolist()
    sample_names = [item for sublist in sample_names for item in sublist]
    sample_names = list(set(sample_names))
    sample_names.sort()

    SNP_seqs = [""] * len(sample_names)

    pos = 0
    for position in df['Position']:
        print(pos)
        # only keep SNPs of impact given in filter (LOW, MODERATE or HIGH)
        sequences = df[df['Position'] == position]["SAMPLES"].str.split(",").tolist()
        sequences = sequences[0]

        for sample in sample_names:
            s = sample_names.index(sample)
            if sample in sequences:
                nuc = df[df['Position'] == position]["ALT_CONTENT"]
                SNP_seqs[s] += nuc
            else:
                nuc = df[df['Position'] == position]["REF_CONTENT"]
                SNP_seqs[s] += nuc
        pos += 1


    print("Mega file of SNP sequences was saved at " + path)
"""

