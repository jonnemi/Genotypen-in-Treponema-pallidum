import math

import pandas as pd
from Bio import SeqIO
import sqlite3
import numpy as np
from collections import Counter
import ast
import os
import dataProcess
import difflib
import csv


def getRefLoci(genome_record):
    # Loop over the features
    data = []
    for feature in genome_record.features:
        if feature.type == "CDS":
            feature_data = {"locus_tag": feature.qualifiers["locus_tag"][0],
                            "location": feature.location,
                            "product": feature.qualifiers["product"][0]}
            data.append(feature_data)
    loci = pd.DataFrame(data)
    return loci


def lociSNVdensity(loci_df):
    # get tSNE input data from snp database
    # establish connection to sqlite database
    conn = sqlite3.connect("snps.db")
    cursor = conn.cursor()

    # get all reference_GenWidePos for SNVs
    sql_command = "SELECT reference_GenWidePos FROM unambiguous ORDER BY reference_GenWidePos;"
    cursor.execute(sql_command)
    content = cursor.fetchall()
    loci_df["var sites"] = loci_df['location'].apply(lambda x: [list(x)[0], list(x)[-1]])
    loci_df["var site count"] = loci_df["var sites"].str.len()
    loci_df["length"] = loci_df["location"].str.len()
    loci_df["SNV density"] = round(loci_df["var site count"] / loci_df["length"], 5)

    sort_df = loci_df.sort_values("SNV density", ascending=False).drop(columns=['location'])
    cols = ["locus_tag", "var site count", "length", "SNV density", "product", 'var sites']
    sort_df = sort_df[cols]
    # save MLST dataframe to file
    sort_df.to_csv("loci.tsv", sep='\t')


def querySNVdensity(loci_df):
    # read tsvFile containing query sequence SNPs
    query = pd.read_csv("variantContentTable.tsv", sep='\t', header=1)
    # create vector of SNPs, with field for each SNP site (reference_GenWidePos)
    sample_names = list(query.columns)
    sample_names.remove('Position')
    sample_names.remove('Reference')

    loci_df["var sites"] = loci_df['location'].map(lambda x: np.intersect1d(x, query["Position"].tolist()))
    loci_df["var site count"] = loci_df["var sites"].str.len()
    loci_df["length"] = loci_df["location"].str.len()
    loci_df["SNV density"] = round(loci_df["var site count"] / loci_df["length"], 5)

    sort_df = loci_df.sort_values("SNV density", ascending=False).drop(columns=['var sites', 'location'])
    cols = ["locus_tag", "var site count", "length", "SNV density", "product"]
    sort_df = sort_df[cols]
    # save MLST dataframe to file
    sort_df.to_csv("loci.tsv", sep='\t')


def rRNAtype(seq_list):
    type_counts = dict.fromkeys(['NA', 'S', 'R8', 'R9'], 0)
    for sequence in seq_list:
        if sequence is None:
            type_counts['NA'] += 1
        else:
            sequence = sequence[0].tolist()
            if 2105 in sequence or 5058 in sequence:
                type_counts['R8'] += 1
            elif 2106 in sequence or 5059 in sequence:
                type_counts['R9'] += 1
            else:
                type_counts['S'] += 1

    total = sum(type_counts.values(), 0.0)
    type_counts = {k: round(v / total, 2) for k, v in type_counts.items()}
    type_counts = {x:y for x,y in type_counts.items() if y!=0.0}

    if len(type_counts) > 1:
        return_string = ""
        for key, value in type_counts.items():
            return_string += key + " (" + str(round(value * 100)) + "%)/"
        return_string = return_string[:-1]
    else:
        return_string = next(iter(type_counts))
    return return_string


def MLST(SNP_vec, ordered_names, all_names, loci_list):
    # dataframe to hold allelic profiles of reference genomes ordered by sequence name
    profiles = pd.DataFrame(all_names, columns=['Name'])
    loci_SNP_cols = []

    # insert locus variant for all loci in loci_list into profiles dataframe
    l = 0
    for locus in SNP_vec:
        uniq = np.unique(locus, axis=0).tolist()
        zipped = list(zip(ordered_names[l], locus))
        col_name = loci_list[l] + "_SNPs"
        loci_SNP_cols.append(col_name)
        df = pd.DataFrame(zipped, columns=['Name', col_name])
        df[loci_list[l]] = df[col_name].apply(lambda x: uniq.index(x.tolist()) + 1)

        profiles = profiles.merge(df, on=['Name'], how='left')

        l += 1
    profiles[loci_list] = profiles[loci_list].astype('Int64')
    del loci_SNP_cols[-1]

    # create new column holding allelic profiles
    subset = loci_list.copy()
    subset.remove("23S rRNA")

    # change NaN entries to X in allelic profile
    profiles[subset] = profiles[subset].astype('str').apply(lambda x: x.replace("<NA>", "X"))
    profiles['Allelic Profile'] = profiles[subset].astype(str).agg('.'.join, axis=1)

    # add column SNP vectors that holds a tuple of all loci SNP vectors defining the allelic profile
    profiles['SNP vectors'] = profiles[loci_SNP_cols].apply(list, axis=1)

    # add column holding the number of samples for each allelic profile
    profiles['No. of samples'] = profiles['Allelic Profile'].map(profiles['Allelic Profile'].value_counts())
    # new column samples holds all samples with given allelic profile instead of one row for each sample
    sample_mapping = profiles.groupby(['Allelic Profile'], as_index=False)['Name'].apply(', '.join)
    profiles['Samples'] = profiles['Allelic Profile'].map(sample_mapping.set_index('Allelic Profile')['Name'])

    # reshape column 23S rRNA_SNPs column to hold all 23S rRNA SNPs for each allelic profile
    profiles['23S rRNA_SNPs'] = profiles['23S rRNA_SNPs'].where(pd.notnull(profiles['23S rRNA_SNPs']), None)
    profiles['23S rRNA_SNPs'] = profiles['23S rRNA_SNPs'].apply(lambda x: x if x is None else np.where(x == "g"))
    rRNA_mapping =  profiles.groupby(['Allelic Profile'])['23S rRNA_SNPs'].apply(list).reset_index(name='23S rRNA_SNPs')
    profiles['23S rRNA_SNPs'] = profiles['Allelic Profile'].map(rRNA_mapping.set_index('Allelic Profile')['23S rRNA_SNPs'])
    profiles['23S rRNA'] = profiles['23S rRNA_SNPs'].apply(lambda x: rRNAtype(x))

    # drop duplicate rows in profiles
    profiles = profiles.drop_duplicates(subset=subset)

    # drop name column
    profiles = profiles.drop('Name', axis=1)
    profiles = profiles.reset_index(drop=True)

    # rearrange columns of profiles dataframe
    cols = profiles.columns.tolist()
    cols.remove('Allelic Profile')
    cols.insert(0, 'Allelic Profile')
    profiles = profiles[cols]

    return profiles


def compareMLST(tsvQuery, loci_list, new_format, filter=[]):
    # get loci training and test data set
    data = dataProcess.getLociDataset(tsvQuery, "snps.db", "string/none", loci_list, new_format, filter)
    print("Dataset acquired")
    print("----------------------------------------------------------------------------")

    # create MLST for referene data set
    ref_profiles = MLST(data['train'], data['train_seq_names'], data['all_train_names'], loci_list)
    ref_profiles['sort'] = ref_profiles['TPANIC_RS00695']
    ref_profiles['sort'] = ref_profiles['sort'].replace("X", "10000")
    ref_profiles['sort'] = pd.to_numeric(ref_profiles['sort'])
    ref_profiles = ref_profiles.sort_values('sort').reset_index(drop=True)
    print("Reference MLST:")
    print(ref_profiles)
    print("Number of allelic profiles: " + str(len(ref_profiles)))
    cols = ["Allelic Profile"]
    cols.extend(loci_list)
    cols.extend(["No. of samples"])
    cols = ['Allelic Profile', 'Samples']
    print_ref_profiles = ref_profiles.copy()
    print_ref_profiles = print_ref_profiles[print_ref_profiles["No. of samples"] > 1]
    print_ref_profiles = print_ref_profiles[cols]
    print_ref_profiles.to_csv('ref_MLST.csv', index=False)

    complete_count = ref_profiles[~ref_profiles['Allelic Profile'].str.contains("X")]
    print("Number of complete allelic profiles: " + str(len(complete_count)))
    print("Total number of samples: " + str(ref_profiles['No. of samples'].sum()))

    s_df = ref_profiles[ref_profiles['23S rRNA'].str.contains("S")]
    na_df = ref_profiles[ref_profiles['23S rRNA'].str.contains("NA")]
    s_count = len(s_df['23S rRNA'].tolist())
    na_count = len(na_df['23S rRNA'].tolist())

    print("Percentage of resistant Allelic Profiles: " + str(s_count + na_count))


    print()
    print("----------------------------------------------------------------------------")
    print()

    # create MLST for query data set
    query_profiles = MLST(data['test'], data['test_seq_names'], data['all_test_names'], loci_list)
    query_profiles['sort1'] = query_profiles['TPANIC_RS00695']
    query_profiles['sort1'] = query_profiles['sort1'].replace("X", "10000")
    query_profiles['sort1'] = pd.to_numeric(query_profiles['sort1'])
    query_profiles['sort2'] = query_profiles['TPANIC_RS02695']
    query_profiles['sort2'] = query_profiles['sort2'].replace("X", "10000")
    query_profiles['sort2'] = pd.to_numeric(query_profiles['sort2'])
    query_profiles = query_profiles.sort_values(['sort1', 'sort2']).reset_index(drop=True)
    print("Query MLST:")
    print(query_profiles)

    print_query_profiles = query_profiles.copy()
    if new_format:
        print_query_profiles = print_query_profiles[print_query_profiles["No. of samples"] > 4]
        #print_query_profiles = print_query_profiles[print_query_profiles["No. of samples"] <= 20]
        print("Allelic profiles representing " + str(round(print_query_profiles["No. of samples"].sum() * 100 /query_profiles["No. of samples"].sum()))
             + "% of samples written to file.")


    print_query_profiles = print_query_profiles[print_query_profiles["No. of samples"] > 1]
    print_query_profiles = print_query_profiles[cols]
    print_query_profiles.to_csv('query_MLST.csv', index=False)

    print("Number of allelic profiles: " + str(len(query_profiles)))
    complete_count = query_profiles[~query_profiles['Allelic Profile'].str.contains("X")]
    print("Number of complete allelic profiles: " + str(len(complete_count)))
    print("Total number of samples: " + str(query_profiles['No. of samples'].sum()))

    s_df = query_profiles[query_profiles['23S rRNA'].str.contains("S")]
    na_df = query_profiles[query_profiles['23S rRNA'].str.contains("NA")]
    s_count = len(s_df['23S rRNA'].tolist())
    na_count = len(na_df['23S rRNA'].tolist())

    print("Percentage of sensitive Allelic Profiles: " + str(s_count + na_count))



    """# compare reference and query allelic profiles
    # turn SNP vectors into tuples of tuples and turn pd series into set
    ref_vecs = set(
        ref_profiles['SNP vectors'].apply(lambda x: tuple([tuple(a) if type(a) is np.ndarray else () for a in x])))
    query_vecs = set(
        query_profiles['SNP vectors'].apply(lambda x: tuple([tuple(a) if type(a) is np.ndarray else () for a in x])))

    # get intersection of both sets
    intersect = ref_vecs.intersection(query_vecs)
    #print(intersect)
    #print(len(intersect))"""

    sample_profiles = pd.DataFrame(data['all_test_names'], columns=['Sample'])
    sample_profiles['Allelic Profile'] = 0
    for index, row in sample_profiles.iterrows():
        sample_profiles['Allelic Profile'][index] = query_profiles[query_profiles['Samples'].str.contains(row['Sample'])]['Allelic Profile'].values[0]

    return(sample_profiles)

    #return sample_profiles


#genome_record = SeqIO.read("NC_021490.2.gb", "genbank")
#lociSNVdensity(getRefLoci(genome_record))
#compareMLST("variantContentTable.tsv", ["TPANIC_RS00695", "TPANIC_RS02695", "TPANIC_RS03500"], False)
#compareMLST("Parr1509_CP004010_SNPSummary.tsv", ["TPANIC_RS00695", "TPANIC_RS02695", "TPANIC_RS03500"], True, ["MODERATE", "HIGH"])