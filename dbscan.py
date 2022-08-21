import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.metrics import v_measure_score
from sklearn.neighbors import NearestNeighbors
from kneed import KneeLocator

def tSNE_DBSCAN(tSNE_tsv):
    # read tSNE data from file
    tSNE = pd.read_csv(tSNE_tsv, sep='\t')
    data = tSNE[['tsne_1', 'tsne_2']].copy()

    # find optimal epsilon for dbscan with elbow method
    nearest_neighbors = NearestNeighbors(n_neighbors=11)
    neighbors = nearest_neighbors.fit(data)

    distances, indices = neighbors.kneighbors(data)
    distances = np.sort(distances[:, 10], axis=0)

    # find elbow/knee with kneed
    i = np.arange(len(distances))
    knee = KneeLocator(i, distances, S=1, curve='convex', direction='increasing', interp_method='polynomial')

    knee.plot_knee()
    plt.xlabel("Points")
    plt.ylabel("Distance")

    epsilon = distances[knee.knee]

    dbscan_cluster1 = DBSCAN(eps=epsilon)
    dbscan_cluster1.fit(data)

    # Visualizing DBSCAN
    fig, ax = plt.subplots(1)
    sns.scatterplot(x=data['tsne_1'], y=data['tsne_2'], hue=dbscan_cluster1.labels_, ax=ax, s=5, style=tSNE['strain'],
                    markers=["o", "v", "D", "X", "P"], palette='Dark2_r')
    ax.set_aspect('equal')
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0, ncol=3)
    fig.tight_layout()

    # Number of Clusters
    labels = dbscan_cluster1.labels_
    N_clus = len(set(labels)) - (1 if -1 in labels else 0)
    print('Estimated no. of clusters: %d' % N_clus)

    # Identify Noise
    n_noise = list(dbscan_cluster1.labels_).count(-1)
    print('Estimated no. of noise points: %d' % n_noise)

    # Calculating v_measure
    print('v_measure =', v_measure_score(tSNE['label'], labels))

    plt.show()


#tSNE_DBSCAN("tSNE.tsv")


def DBSCAN_cluster(df):
    data = df.iloc[:, 0:2]
    min_samples = 10

    # find optimal epsilon for dbscan with elbow method
    nearest_neighbors = NearestNeighbors(n_neighbors=min_samples)
    neighbors = nearest_neighbors.fit(data)

    distances, indices = neighbors.kneighbors(data)
    distances = np.sort(distances[:, min_samples-1], axis=0)

    # find elbow/knee with kneed
    i = np.arange(len(distances))
    knee = KneeLocator(i, distances, S=1, curve='convex', direction='increasing', interp_method='polynomial')

    #knee.plot_knee()
    #plt.xlabel("Points")
    #plt.ylabel("Distance")

    epsilon = distances[knee.knee]

    dbscan_cluster1 = DBSCAN(eps=0.9, min_samples=3)
    #dbscan_cluster1 = DBSCAN(eps=epsilon, min_samples=min_samples)
    dbscan_cluster1.fit(data)

    # Number of Clusters
    labels = dbscan_cluster1.labels_
    N_clus = len(set(labels)) - (1 if -1 in labels else 0)

    df['cluster'] = labels

    cluster_centers = pd.DataFrame(set(labels), columns=['cluster'])

    cluster_seqs = df.groupby('cluster', as_index=False)['label'].apply(', '.join)
    cluster_centers['cluster_sequences'] = cluster_centers['cluster'].map(
        cluster_seqs.set_index('cluster')['label'])
    cluster_centers['cluster_sequences'] = cluster_centers['cluster_sequences'].str.split(",")
    cluster_centers['size'] = cluster_centers['cluster_sequences'].apply(lambda x: len(x))

    return cluster_centers
