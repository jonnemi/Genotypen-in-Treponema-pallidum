import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import scipy.cluster.hierarchy as shc
from sklearn.cluster import AgglomerativeClustering

# read tSNE data from file
tSNE = pd.read_csv("tSNE.tsv", sep='\t')
data = tSNE[['tsne_1', 'tsne_2']].copy()

# create dendrogram for dataset using scipy
plt.figure(figsize=(10, 7))
plt.title("Dendogram")
dend = shc.dendrogram(shc.linkage(data, method='ward'))

# cluster data with agglomerative clustering
cluster = AgglomerativeClustering(n_clusters=5, affinity='euclidean', linkage='ward')
cluster.fit_predict(data)

# plot hierarchical clustering of tSNE
fig, ax = plt.subplots(1)
sns.scatterplot(x=data['tsne_1'], y=data['tsne_2'], hue=cluster.labels_, ax=ax, s=5, style=tSNE['strain'],
                markers=["o", "v", "D", "X", "P"], palette='Dark2_r')
ax.set_aspect('equal')
ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0, ncol=3)
fig.tight_layout()
plt.show()