# Biopython's SeqIO module handles sequence input/output
import pandas as pd
from Bio import SeqIO
import sqlite3
import numpy as np
from collections import Counter


def getRefLoci(seq_record):
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
    loci_df["var site count"] = loci_df["var sites"].str.len()
    loci_df["length"] = loci_df["location"].str.len()
    loci_df["SNV density"] = loci_df["var site count"] / loci_df["length"]

    sort_df = loci_df.sort_values('var site count', ascending=False)

    # save MLST dataframes to file
    sort_df.to_csv("MLST.tsv", sep='\t')

genome_record = SeqIO.read("NC_021490.2.gb", "genbank")
lociSNVdensity(getRefLoci(genome_record))