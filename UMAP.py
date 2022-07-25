import numpy as np
from sklearn.datasets import load_digits
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import dataProcess
import umap


def snpUmap(queryTSV, db_name):
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

snpUmap("variantContentTable.tsv", "snps.db")
