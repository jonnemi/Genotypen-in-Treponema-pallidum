import pandas as pd
from Bio import SeqIO
import sqlite3
import numpy as np
from collections import Counter
import ast
import os
import dataProcess
import difflib


def getRefLoci(genome_record):
    # Loop over the features
    data = []
    for feature in genome_record.features:
        if feature.type == "CDS":
            feature_data = {"locus_tag": feature.qualifiers["locus_tag"][0],
                            "location": list(feature.location),
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
    sites = list()
    for pos in content:
        sites.append(pos[0])
    conn.close()

    loci_df["var sites"] = loci_df['location'].map(lambda x: np.intersect1d(x, sites))
    site_count = Counter(sites)
    loci_df["var site count"] = loci_df["var sites"].str.len()
    loci_df["length"] = loci_df["location"].str.len()
    loci_df["SNV density"] = round(loci_df["var site count"] / loci_df["length"], 5)

    sort_df = loci_df.sort_values("SNV density", ascending=False).drop(columns=['location'])
    cols = ["locus_tag", "var site count", "length", "SNV density", "product", 'var sites']
    sort_df = sort_df[cols]
    sort_df["var sites"] = sort_df["var sites"].apply(lambda x: x.tolist())
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


def MLST(SNP_vec, ordered_names, all_names, loci_list):
    # dataframe to hold allelic profiles of reference genomes ordered by sequence name
    profiles = pd.DataFrame(all_names, columns=['Name'])
    # insert locus variant for all loci in loci_list into profiles dataframe
    l = 0
    for loci in SNP_vec:
        uniq = np.unique(loci, axis=0).tolist()
        zipped = list(zip(ordered_names[l], loci))
        df = pd.DataFrame(zipped, columns=['Name', 'SNPs'])
        df[loci_list[l]] = df['SNPs'].apply(lambda x: uniq.index(x.tolist()) + 1)

        profiles = profiles.merge(df, on=['Name'], how='left')
        profiles = profiles.drop('SNPs', axis=1)

        l += 1
    profiles[loci_list] = profiles[loci_list].astype('Int64')


    # create new column holding allelic profiles
    subset = loci_list.copy()
    subset.remove("23S rRNA")
    profiles['Allelic Profile'] = profiles[subset].astype(str).agg('.'.join, axis=1)

    # change NaN entries to X in allelic profile
    profiles['Allelic Profile'] = profiles['Allelic Profile'].astype('str').apply(lambda  x: x.replace("<NA>", "X"))

    # add column holding the number of samples for each allelic profile
    profiles['No. of samples'] = profiles['Allelic Profile'].map(profiles['Allelic Profile'].value_counts())

    # new column samples holds all samples with given allelic profile instead of one row for each sample
    profiles['samples'] = profiles['Allelic Profile'].apply(lambda x: profiles[profiles['Allelic Profile']==x]['Name'].tolist())

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

    print(profiles)



def compareMLST(tsvQuery, loci_list):
    # get loci training and test data set
    data = dataProcess.getLociDataset(tsvQuery, "snps.db", ".", loci_list)

    # create MLST for referene data set
    ref_profiles = MLST(data['train'], data['train_seq_names'], data['all_train_names'], loci_list)

    # create MLST for query data set
    query_profiles = MLST(data['test'], data['test_sample_names'], data['test_sample_names'], loci_list)

# genome_record = SeqIO.read("NC_021490.2.gb", "genbank")
# lociSNVdensity(getRefLoci(genome_record))
compareMLST("variantContentTable.tsv", ["TPANIC_RS00695", "TPANIC_RS02695", "TPANIC_RS03500"])