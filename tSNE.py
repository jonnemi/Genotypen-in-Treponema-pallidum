import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os

from openTSNE import TSNE
import dataProcess


def tSNE(queryTSV, db_name):
    # load train and test dataset
    dataset = dataProcess.getADAdataset(queryTSV, db_name, [0, 0, 0, 0, 1])

    # fit and transform standard tSNE
    fit_standard = TSNE(
        perplexity=30,
        initialization="random",
        metric="euclidean",
        n_jobs=8,
        random_state=3,
    ).fit(dataset["train"])
    embedding_standard = fit_standard.transform(dataset["test"])

    # fit and transform tSNE with PCA initialization
    fit_pca = TSNE(
        perplexity=30,
        initialization="pca",
        metric="euclidean",
        n_jobs=8,
        random_state=3,
    ).fit(dataset["train"])
    embedding_pca = fit_pca.transform(dataset["test"])

    # fit and transform tSNE using cosine distance
    fit_cosine = TSNE(
        perplexity=30,
        initialization="random",
        metric="hamming",
        n_jobs=8,
        random_state=3,
    ).fit(dataset["train"])
    embedding_cosine = fit_cosine.transform(dataset["test"])

    # fit and transform tSNE with PCA initialization using cosine distance
    fit_pca_cosine = TSNE(
        perplexity=30,
        initialization="pca",
        metric="hamming",
        n_jobs=8,
        random_state=3,
    ).fit(dataset["train"])
    embedding_pca_cosine = fit_pca_cosine.transform(dataset["test"])

    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(17, 17))

    titles = ["Standard t-SNE", "PCA initialization", "Hamming distance", "PCA initialization + Hamming distance"]
    fits = [fit_standard, fit_pca, fit_cosine, fit_pca_cosine]
    embeddings = [embedding_standard, embedding_pca, embedding_cosine, embedding_pca_cosine]

    strains = dataProcess.assignStrains(dataset['test_sample_names'])

    j = [0, 0, 1, 1]
    for i in range(0, 4):
        # plot hierarchical clustering of tSNE

        train_pca_df = pd.DataFrame(
            {'tsne_1': fits[i][:, 0], 'tsne_2': fits[i][:, 1],
             'label': dataset['train_seq_names']})
        test_pca_df = pd.DataFrame(
            {'tsne_1': embeddings[i][:, 0], 'tsne_2': embeddings[i][:, 1],
             'label': dataset['test_sample_names']})
        test_pca_df["strain"] = strains

        n = i % 2
        m = j[i]
        sns.scatterplot(x='tsne_1', y='tsne_2', data=train_pca_df, ax=axs[n, m], s=5, alpha=0.5,
                        color='black')
        sns.scatterplot(x='tsne_1', y='tsne_2', hue='strain', data=test_pca_df, ax=axs[n, m], s=5)
        axs[n, m].set_aspect('equal')
        axs[n, m].set_xticklabels([])
        axs[n, m].set_yticklabels([])
        # axs[n, m].set_xlabel('')
        # axs[n, m].set_ylabel('')
        axs[n, m].set_title(titles[i])
        axs[n, m].legend([], [], frameon=False)
    plt.gca().get_legend_handles_labels()
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)
    plt.show()


def tuneTSNE(queryTSV, db_name):
    # load train and test dataset
    dataset = dataProcess.getADAdataset(queryTSV, db_name, [0, 0, 0, 0, 1])
    strains = dataProcess.assignStrains(dataset['test_sample_names'])
    fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(17, 17))
    perplexities = [2, 5, 10, 20, 30, 50]
    j = [0, 0, 0, 1, 1, 1]

    # fit and transform tSNE with PCA initialization using cosine distance with different perplexities
    for i in range(0, 6):
        fit = TSNE(
            perplexity=perplexities[i],
            initialization="pca",
            metric="hamming",
            n_jobs=8,
            random_state=3,
        ).fit(dataset["train"])
        embedding = fit.transform(dataset["test"])

        title = "Perplexity = " + str(perplexities[i])

        # plot hierarchical clustering of tSNE
        train_pca_df = pd.DataFrame(
            {'tsne_1': fit[:, 0], 'tsne_2': fit[:, 1],
             'label': dataset['train_seq_names']})
        test_pca_df = pd.DataFrame(
            {'tsne_1': embedding[:, 0], 'tsne_2': embedding[:, 1],
             'label': dataset['test_sample_names']})
        test_pca_df["strain"] = strains

        n = i % 3
        m = j[i]
        sns.scatterplot(x='tsne_1', y='tsne_2', data=train_pca_df, ax=axs[n, m], s=5, alpha=0.5,
                        color='black')
        sns.scatterplot(x='tsne_1', y='tsne_2', hue='strain', data=test_pca_df, ax=axs[n, m], s=5)
        axs[n, m].set_aspect('equal')
        axs[n, m].set_xticklabels([])
        axs[n, m].set_yticklabels([])
        axs[n, m].set_title(title)
        axs[n, m].legend([], [], frameon=False)
    plt.gca().get_legend_handles_labels()
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)
    plt.show()


def compareTrainTest(queryTSV, db_name):
    # load train and test dataset
    dataset = dataProcess.getADAdataset(queryTSV, db_name, [0, 0, 0, 0, 1])

    fig, ax = plt.subplots(1)

    # fit and transform tSNE with PCA initialization using cosine distance with different perplexities
    fit = TSNE(
        perplexity=30,
        initialization="pca",
        metric="hamming",
        n_jobs=8,
        random_state=3,
    ).fit(dataset["train"])
    embedding = fit.transform(dataset["test"])

    # plot hierarchical clustering of tSNE
    train_pca_df = pd.DataFrame(
        {'tsne_1': fit[:, 0], 'tsne_2': fit[:, 1],
         'label': dataset['train_seq_names']})
    test_pca_df = pd.DataFrame(
        {'tsne_1': embedding[:, 0], 'tsne_2': embedding[:, 1],
         'label': dataset['test_sample_names']})

    names = dataProcess.getStrains(dataset['train_seq_names'])
    train_pca_df['label'] = names
    intersection = list(set(dataset['test_sample_names']).intersection(names))
    print(len(intersection))

    train_pca_df = train_pca_df[train_pca_df['label'].isin(test_pca_df['label'])].reset_index()
    test_pca_df = test_pca_df[test_pca_df['label'].isin(intersection)].reset_index()

    sns.scatterplot(x='tsne_1', y='tsne_2', data=train_pca_df, ax=ax, s=5, alpha=0.5,
                    hue='label')
    sns.scatterplot(x='tsne_1', y='tsne_2', hue='label', data=test_pca_df, ax=ax, s=5)
    ax.set_aspect('equal')
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.legend([], [], frameon=False)
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0, ncol=3)
    plt.show()


#tSNE("variantContentTable.tsv", "snps.db")
#tuneTSNE("variantContentTable.tsv", "snps.db")
#compareTrainTest("variantContentTable.tsv", "snps.db")



def querytSNE(tsvFile, filter, default_enc):
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

    # create encoded SNP vector, with field for each SNP site using one-hot encoding
    process = dataProcess.newEncQuery(query, default_enc, positions)
    enc_2D = process["data"]
    sample_names = process["sample_names"]

    fit = TSNE(
        perplexity=30,
        initialization="pca",
        metric="hamming",
        n_jobs=8,
        random_state=3,
        verbose=True,
    ).fit(enc_2D)

    tsne_df = pd.DataFrame(
        {'tsne_1': fit[:, 0], 'tsne_2': fit[:, 1],
         'label': sample_names})

    fig, ax = plt.subplots(1)
    sns.scatterplot(x='tsne_1', y='tsne_2', data=tsne_df, ax=ax, s=5, hue='label')
    ax.set_aspect('equal')
    ax.legend([],[], frameon=False)
    plt.show()

#querytSNE("Parr1509_CP004010_SNPSummary.tsv", ["MODERATE", "HIGH"], "one-hot")