import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import seaborn as sns
import numpy as np
from kneed import KneeLocator

# read tSNE data from file
tSNE = pd.read_csv("tSNE.tsv", sep='\t')
data = tSNE[['tsne_1', 'tsne_2']].copy()

# run k-Means with a range of k
distortions = []
K = range(1, 20)
for k in K:
    kmeanModel = KMeans(n_clusters=k)
    kmeanModel.fit(data)
    distortions.append(kmeanModel.inertia_)

# plot distortions of k-Means
plt.figure(figsize=(16, 8))
plt.plot(K, distortions, 'bx-')
plt.xlabel('k')
plt.ylabel('Distortion')
plt.title('The Elbow Method showing the optimal k')
plt.xticks(K)

# elbow can be observed for k=3, fit model for k=3
kmeanModel = KMeans(n_clusters=3)
kmeanModel.fit(data)

tSNE['k-means'] = kmeanModel.predict(data)

# plot k-Means clustering of tSNE
fig, ax = plt.subplots(1)
sns.scatterplot(x='tsne_1', y='tsne_2', hue='k-means', data=tSNE, ax=ax, s=5, style='strain',
                markers=["o", "v", "D", "X", "P"], palette='Dark2_r')
ax.set_aspect('equal')
ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0, ncol=3)
fig.tight_layout()
plt.show()
