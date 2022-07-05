from sklearn.decomposition import PCA
import pandas as pd
import matplotlib.pyplot as plt
import dataProcess
import seaborn as sns


def snpPCA(queryTSV, db_name):
    dataset = dataProcess.getADAdataset(queryTSV, db_name, [0, 0, 0, 0, 1, 0])
    pca = PCA(n_components=2)
    pca_fit = pca.fit(dataset["train"])
    trainPC = pca.transform(dataset["train"])
    testPC = pca.transform(dataset["test"])
    strains = dataProcess.assignStrains(dataset['test_sample_names'])

    train_pca_df = pd.DataFrame(
        {'principal_comp_1': trainPC[:, 0], 'principal_comp_2': trainPC[:, 1],
         'label': dataset['train_seq_names']})
    test_pca_df = pd.DataFrame(
        {'principal_comp_1': testPC[:, 0], 'principal_comp_2': testPC[:, 1],
         'label': dataset['test_sample_names']})
    test_pca_df["strain"] = strains
    fig, ax = plt.subplots(1)
    sns.scatterplot(x='principal_comp_1', y='principal_comp_2', data=train_pca_df, ax=ax, s=5, alpha=0.5, color='black')
    sns.scatterplot(x='principal_comp_1', y='principal_comp_2', hue='label', data=test_pca_df, ax=ax, s=5, style='strain', markers=["o", "v", "D", "X"])
    ax.set_aspect('equal')
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0, ncol=3)
    fig.tight_layout()
    plt.show()


snpPCA("variantContentTable.tsv", "snps.db")
