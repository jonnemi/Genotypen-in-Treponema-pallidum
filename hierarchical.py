import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import scipy.cluster.hierarchy as shc
from sklearn.cluster import AgglomerativeClustering

# read tSNE data from file
tSNE = pd.read_csv("tSNE.tsv", sep='\t')
data = tSNE[['tsne_1', 'tsne_2']].copy()

# create dendrogram for dataset using scipy
"""plt.figure(figsize=(10, 7))
plt.title("Dendogram")
dend = shc.dendrogram(shc.linkage(data, method='single'))"""

# cluster data with agglomerative clustering
linkage_list = ("ward", "complete", "average", "single")
fig, axs = plt.subplots(2, 2)

for i in range(0, 3):
    cluster = AgglomerativeClustering(n_clusters=5, affinity='euclidean', linkage=linkage_list[i])
    cluster.fit_predict(data)
    title = linkage_list[i] + " linkage"

    if i == 0:
        j = 0
    elif i == 1:
        j = 0
    elif i == 2:
        j = 1
        i = 0
    else:
        j = 1
        i = 1

    # plot hierarchical clustering of tSNE

    sns.scatterplot(x=data['tsne_1'], y=data['tsne_2'], hue=cluster.labels_, ax=axs[i, j], s=5, palette='Dark2_r')
    axs[i, j].set_aspect('equal')
    axs[i, j].set_title(title)
    axs[i, j].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0, ncol=3)
    fig.tight_layout()
plt.show()
