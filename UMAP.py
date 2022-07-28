import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import dataProcess
import umap
import os


def adaUmap(queryTSV, db_name):
    dataset = dataProcess.getADAdataset(queryTSV, db_name, [0, 0, 0, 0, 1])

    reducer = umap.UMAP()
    reducer.fit(dataset["train"])
    train= reducer.transform(dataset["train"])
    test = reducer.transform(dataset["test"])

    strains = dataProcess.assignStrains(dataset['test_sample_names'])

    train_pca_df = pd.DataFrame(
        {'umap_1': train[:, 0], 'umap_2': train[:, 1],
         'label': dataset['train_seq_names']})
    test_pca_df = pd.DataFrame(
        {'umap_1': test[:, 0], 'umap_2': test[:, 1],
         'label': dataset['test_sample_names']})
    test_pca_df["strain"] = strains
    fig, ax = plt.subplots(1)
    sns.scatterplot(x='umap_1', y='umap_2', data=train_pca_df, ax=ax, s=5, alpha=0.5, color='black')
    sns.scatterplot(x='umap_1', y='umap_2', hue='strain', data=test_pca_df, ax=ax, s=5, style='strain', markers=["D", "v"])
    ax.set_aspect('equal')
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0, ncol=3)
    fig.tight_layout()
    plt.show()

#adaUmap("variantContentTable.tsv", "snps.db")


def queryUmap(tsvFile, filter, default_enc):
    if os.path.isfile(tsvFile):
        query = pd.read_csv(tsvFile, sep='\t', header=0)
    query.rename(columns={'POSITION': 'Position'}, inplace=True)

    # create vector of SNPs, with field for each SNP site
    # remove all indel rows from dataframe
    query = query[query['ALT_CONTENT'].str.len() == 1 & ~query['ALT_CONTENT'].str.contains("-")]
    query.reset_index()

    # get all positions as list
    positions = query["Position"].tolist()

    # create encoded SNP vector, with field for each SNP site using one-hot encoding
    process = dataProcess.newEncQuery(query, default_enc, positions, filter)
    enc_2D = process["data"]
    sample_names = process["sample_names"]


    reducer = umap.UMAP(metric='hamming')
    transform = reducer.fit_transform(enc_2D)

    umap_df = pd.DataFrame(
        {'umap_1': transform[:, 0], 'umap_2': transform[:, 1],
         'label': sample_names})

    fig, ax = plt.subplots(1)
    sns.scatterplot(x='umap_1', y='umap_2', data=umap_df, ax=ax, s=5, hue='label')
    ax.set_aspect('equal')
    ax.legend([],[], frameon=False)
    plt.show()

queryUmap("Parr1509_CP004010_SNPSummary.tsv", ["MODERATE", "HIGH"], [0, 0, 0, 0, 1])
