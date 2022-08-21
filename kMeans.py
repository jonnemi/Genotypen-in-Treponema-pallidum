import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import seaborn as sns
import numpy as np
from kneed import KneeLocator

def tSNE_KMeans(tSNE_tsv):
    # read tSNE data from file
    tSNE = pd.read_csv(tSNE_tsv, sep='\t')
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

#tSNE_KMeans("tSNE.tsv")

def KMeans_cluster(df):
    data = df.iloc[:, 0:2]

    # run k-Means with a range of k
    distortions = []
    K = range(1, 100)
    for k in K:
        kmeanModel = KMeans(n_clusters=k)
        kmeanModel.fit(data)
        distortions.append(kmeanModel.inertia_)

    # elbow can be observed for k=3, fit model for k=3
    knee = KneeLocator(K, distortions, S=1, curve='convex', direction='decreasing', interp_method='polynomial')
    k = K[knee.knee]
    #knee.plot_knee()

    kmeanModel = KMeans(n_clusters=k)
    kmeanModel.fit(data)

    df['cluster'] = kmeanModel.predict(data)

    # create dataframe with with cluster number and cluster center location
    cluster_centers = pd.DataFrame(kmeanModel.cluster_centers_, columns=['UMAP_1', 'UMAP_2'])
    cluster_centers['cluster'] = list(range(0, k))

    cluster_seqs = df.groupby('cluster', as_index=False)['label'].apply(', '.join)
    cluster_centers['cluster_sequences'] = cluster_centers['cluster'].map(cluster_seqs.set_index('cluster')['label'])
    cluster_centers['cluster_sequences'] = cluster_centers['cluster_sequences'].str.split(",")
    cluster_centers['size'] = cluster_centers['cluster_sequences'].apply(lambda x: len(x))

    return cluster_centers

