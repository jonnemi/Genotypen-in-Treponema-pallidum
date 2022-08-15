import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import dataProcess
import umap
import os
import kMeans
from openTSNE import TSNE
import ast



def compare(tsvFile, filter, default_enc, loci=list()):
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

    # get all positions as list
    positions = query["Position"].tolist()

    # check if MLST is wanted
    if loci:
        # load loci dataframe
        loci_df = pd.read_csv("loci.tsv", sep='\t', converters={"var sites": ast.literal_eval})
        # create list of MLST loci SNP positions
        MLST_loci = loci_df[loci_df['locus_tag'].isin(loci)]
        MLST_positions = []
        MLST_loci["var sites"].apply(lambda x: MLST_positions.extend(x))

        rRNA_R8_v1 = 235204
        rRNA_R9_v1 = 235205
        rRNA_R8_v2 = 283649
        rRNA_R9_v2 = 283650

        rRNA_pos = [rRNA_R8_v1, rRNA_R9_v1, rRNA_R8_v2, rRNA_R9_v2]

        MLST_positions.extend(rRNA_pos)
        MLST_positions.sort()

        query = query[query["Position"].isin(MLST_positions)].reset_index()

    # create encoded SNP vector, with field for each SNP site using the given encoding
    process = dataProcess.newEncQuery(query, default_enc, positions)
    enc_2D = process["data"]
    sample_names = process["sample_names"]

    print("Dataset acquired, starting UMAP...")

    # apply UMAP to data
    reducer = umap.UMAP(metric='manhattan',
                        n_neighbors = 30,
                        min_dist = 0.0,
                        n_components = 2,
                        random_state = 42,
                        verbose=True
                        )
    transform = reducer.fit_transform(enc_2D)

    # process umap dimension reduction with kMeans
    umap_df = pd.DataFrame(
        {'UMAP_1': transform[:, 0], 'UMAP_2': transform[:, 1],
         'label': sample_names})


    # define clusters in UMAP projection using k-means
    k_Means = kMeans.UMAP_KMeans(umap_df)

    # dataframe holding kMeans cluster label for each sequence
    sample_clustering = pd.DataFrame(sample_names, columns=['sample_names'])
    sample_clustering['cluster'] = -1
    c = 0
    for cluster in k_Means['cluster_sequences']:
        cluster = [x.replace(" ", "") for x in cluster]
        sample_clustering.loc[sample_clustering['sample_names'].isin(cluster), ['cluster']] = c
        c += 1



    fig, ax = plt.subplots(1, 2)
    sns.scatterplot(x='UMAP_1', y='UMAP_2', data=umap_df, ax=ax[0], s=5, hue='K-means_cluster', palette='colorblind')
    ax[0].legend([], [], frameon=False)

    # fit TSNE as compairson to UMAP
    fit = TSNE(
        perplexity=30,
        initialization="pca",
        metric="manhattan",
        n_jobs=8,
        random_state=3,
        verbose=True,
    ).fit(enc_2D)

    tsne_df = umap_df.copy()
    tsne_df['tSNE_1'] = fit[:, 0]
    tsne_df['tSNE_2'] = fit[:, 1]

    sns.scatterplot(x='tSNE_1', y='tSNE_2', data=tsne_df, ax=ax[1], s=5, hue='K-means_cluster', palette='colorblind')
    ax[1].legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, title='K-means cluster')
    plt.subplots_adjust(wspace=0.5)
    plt.tight_layout()


    plt.show()


compare("Parr1509_CP004010_SNPSummary.tsv", ["MODERATE", "HIGH"], "binary")