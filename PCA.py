from sklearn.decomposition import PCA
import pandas as pd
import matplotlib.pyplot as plt
import dataProcess
import seaborn as sns


def snpPCA(queryTSV, db_name):
    dataset = dataProcess.getADAdataset(queryTSV, db_name, [0, 0, 0, 0, 1, 0])
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(dataset["train"])

    train_pca_df = pd.DataFrame(
        {'principal_comp_1': principalComponents[:, 0], 'principal_comp_2': principalComponents[:, 1],
         'label': dataset['train_seq_names']})
    fig, ax = plt.subplots(1)
    sns.scatterplot(x='principal_comp_1', y='principal_comp_2', data=train_pca_df, ax=ax, s=5, alpha=0.5, color="black")
    ax.set_aspect('equal')
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0, ncol=3)
    fig.tight_layout()
    plt.show()


snpPCA("variantContentTable.tsv", "snps.db")
