import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import dataProcess
import umap
import os
import kMeans
from openTSNE import TSNE
import ast
import dbscan
import MLST
from sklearn.metrics.cluster import adjusted_rand_score


def legend_cols(cluster_col):
    cluster_count = cluster_col.nunique()
    cols = int(round(cluster_count / 25))
    if cols == 0:
        cols = 1
    return cols


def plotCompare(df, title, umap_titel, tsne_titel):
    fig, ax = plt.subplots(2, 1)
    fig.set_size_inches(11.69,8.27)
    fig.subplots_adjust(top=0.8)
    sns.scatterplot(x='UMAP_1', y='UMAP_2', data=df, ax=ax[0], s=15, hue='umap_cluster',
                    palette=sns.color_palette("hls", df['umap_cluster'].nunique()))
    ax[0].set_title(umap_titel)
    # ax[0].legend([], [], frameon=False)
    ax[0].legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, title='UMAP Cluster',
                 ncol=legend_cols(df['umap_cluster']))

    sns.scatterplot(x='tSNE_1', y='tSNE_2', data=df, ax=ax[1], s=15, hue='tsne_cluster',
                    palette=sns.color_palette("hls", df['tsne_cluster'].nunique()))
    ax[1].legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, title='tSNE Cluster',
                 ncol=legend_cols(df['tsne_cluster']))
    # ax[1].legend([], [], frameon=False)
    ax[1].set_title(tsne_titel)
    plt.subplots_adjust(wspace=0.5)
    plt.tight_layout()
    fig.suptitle(title, y=0.98)


def plotCompareEmbedding(df1, rep_umap, rep_tsne, title):
    """fig, ax = plt.subplots(1, 2)
    fig.subplots_adjust(top=0.8)
    sns.scatterplot(x='UMAP_1', y='UMAP_2', data=df1, ax=ax[0], s=5, hue='umap_cluster',
                    palette=sns.color_palette("hls", df1['umap_cluster'].nunique()))
    sns.scatterplot(x='UMAP_1', y='UMAP_2', data=rep_umap, ax=ax[0], s=10,
                    palette=sns.color_palette("hls", df1['umap_cluster'].nunique()), markers=['D'])
    ax[0].set_title("UMAP")
    # ax[0].legend([], [], frameon=False)
    ax[0].legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, title='UMAP Cluster',
                 ncol=legend_cols(df1['umap_cluster']))

    sns.scatterplot(x='tSNE_1', y='tSNE_2', data=df1, ax=ax[1], s=5, hue='tsne_cluster',
                    palette=sns.color_palette("hls", df1['tsne_cluster'].nunique()))
    sns.scatterplot(x='tSNE_1', y='tSNE_2', data=rep_tsne, ax=ax[1], s=10,
                    palette=sns.color_palette("hls", df1['tsne_cluster'].nunique()), markers=['D'])
    ax[1].legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, title='tSNE Cluster',
                 ncol=legend_cols(df1['tsne_cluster']))
    # ax[1].legend([], [], frameon=False)
    ax[1].set_title("tSNE")
    plt.subplots_adjust(wspace=0.5)
    plt.tight_layout()
    fig.suptitle(title, y=0.98)"""
    rep_umap['Style'] = "Cluster Vektor"

    fig, ax = plt.subplots(1)
    fig.subplots_adjust(top=0.8)
    sns.scatterplot(x='UMAP_1', y='UMAP_2', data=df1, ax=ax, s=5, hue='umap_cluster', hue_order=df1['umap_cluster'].unique().tolist().sort(),
                    palette=sns.color_palette("hls", df1['umap_cluster'].nunique()))
    sns.scatterplot(x='UMAP_1', y='UMAP_2', data=rep_umap, ax=ax, s=15, color='black', markers=['X'], style='Style')
    ax.set_title(title)
    # ax[0].legend([], [], frameon=False)
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, title='UMAP Cluster',
              ncol=legend_cols(df1['umap_cluster']))


def compare(tsvFile, filter, default_enc, clustering, loci=list()):
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

    # check if MLST is wanted
    if loci:
        # load loci dataframe
        loci_df = pd.read_csv("loci.tsv", sep='\t', converters={"var sites": ast.literal_eval})
        # create list of MLST loci SNP positions
        MLST_loci = loci_df[loci_df['locus_tag'].isin(loci)]
        MLST_positions = []
        MLST_loci["var sites"].apply(lambda x: MLST_positions.extend(list(range(min(x) + 1, max(x) + 1))))
        MLST_positions.sort()

        query = query[query["Position"].isin(MLST_positions)].reset_index()

    # create encoded SNP vector, with field for each SNP site using the given encoding
    process = dataProcess.newEncQuery(query, default_enc, positions)
    enc_2D = process["data"]
    sample_names = process["sample_names"]

    print("Dataset acquired, starting UMAP...")

    # apply UMAP to data
    reducer = umap.UMAP(metric='manhattan',
                        n_neighbors=30,
                        min_dist=0.0,
                        n_components=2,
                        random_state=42,
                        verbose=True
                        )
    transform = reducer.fit_transform(enc_2D)

    umap_df = pd.DataFrame(
        {'UMAP_1': transform[:, 0], 'UMAP_2': transform[:, 1],
         'label': sample_names})

    # fit TSNE as compairson to UMAP
    fit = TSNE(
        perplexity=30,
        initialization="pca",
        metric="manhattan",
        n_jobs=8,
        random_state=3,
        verbose=True,
    ).fit(enc_2D)

    tsne_df = pd.DataFrame(
        {'tSNE_1': fit[:, 0], 'tSNE_2': fit[:, 1],
         'label': sample_names})

    # define clusters in UMAP and tSNE projection using given clustering method
    print()
    print("Compute clustering in projections using " + clustering + "...")
    if clustering == "kMeans":
        umap_clusters = kMeans.KMeans_cluster(umap_df)
        print("UMAP clustering done.")
        tsne_clusters = kMeans.KMeans_cluster(tsne_df)
        print("tSNE clustering done.")
    elif clustering == "dbscan":
        umap_clusters = dbscan.DBSCAN_cluster(umap_df)
        print("UMAP clustering done.")
        tsne_clusters = dbscan.DBSCAN_cluster(tsne_df)
        print("tSNE clustering done.")
    else:
        umap_df['cluster'] = -1
        tsne_df['cluster'] = -1
    print("Clustering complete.")

    # define genotypes by SNP frequencies in each cluster
    print()
    print("Compute UMAP genotype by SNP-frequencies...")
    umap_cluster_represent = freqGenotype(umap_df, enc_2D, "umap_SNP_freqs.csv")
    print()
    print("Compute tSNE genotype by SNP-frequencies...")
    tsne_cluster_represent = freqGenotype(tsne_df, enc_2D, "tsne_SNP_freqs.csv")
    print()

    rep_transform = reducer.transform(umap_cluster_represent)

    rep_umap = pd.DataFrame(
        {'UMAP_1': rep_transform[:, 0], 'UMAP_2': rep_transform[:, 1],
         'label': umap_df['cluster'].unique().tolist()})

    rep_fit = fit.transform(tsne_cluster_represent)
    rep_tsne = pd.DataFrame(
        {'tSNE_1': rep_fit[:, 0], 'tSNE_2': rep_fit[:, 1],
         'label': tsne_df['cluster'].unique().tolist()})

    # combine UMAP and tSNE df into one df
    umap_df.rename({'cluster': 'umap_cluster'}, axis=1, inplace=True)
    tsne_df.rename({'cluster': 'tsne_cluster'}, axis=1, inplace=True)
    ML_df = pd.merge(umap_df, tsne_df, how='outer', on='label')
    title = "UMAP und tSNE mit " + clustering + "-Clustering"
    plotCompareEmbedding(ML_df, rep_umap, rep_tsne, title)

    return ML_df


def clusterCenterGenotype(umap_df, k_Means, reducer, enc_2D):
    # inverse transform kMeans cluster centers into high dimensional space
    cluster_centers = k_Means[['umap_1', 'umap_2']].to_numpy()
    inv_transformed_points = reducer.inverse_transform(cluster_centers)

    # tune threshold that bins continuous features of high dimensional cluster centers into binary
    tune_threshold = [0.001, 0.001, 0.01]
    tune_threshold.extend(np.linspace(0, 1, 11))

    errors = np.empty(len(tune_threshold))

    for threshold in tune_threshold:
        # round cluster centers vector entries to either 0 or 1
        p = 0
        inv_transformed_centers = inv_transformed_points.copy()
        for point in inv_transformed_centers:
            point = np.where(point > threshold, 1, 0)
            inv_transformed_centers[p] = point
            p += 1

        # determine defining SNPs of each cluster center vector
        p = 0
        def_positions = []
        for point in inv_transformed_centers:
            other_points = np.delete(inv_transformed_centers, p, axis=0)
            for other in other_points:
                comp = np.equal(point, other)
                diff_pos = np.where(comp == False)[0].tolist()
                def_positions.extend(diff_pos)
            p += 1

        def_positions = list(set(def_positions))
        def_positions.sort()

        def_cluster_snps = []
        for point in inv_transformed_centers:
            point_snps = [point[i] for i in def_positions]
            def_cluster_snps.append(point_snps)

        # dataframe holding kMeans cluster label for each sequence
        sample_names = umap_df['label'].tolist()
        sample_clustering = pd.DataFrame(sample_names, columns=['sample_names'])
        sample_clustering['cluster'] = -1
        c = 0
        for cluster in k_Means['cluster_sequences']:
            cluster = [x.replace(" ", "") for x in cluster]
            sample_clustering.loc[sample_clustering['sample_names'].isin(cluster), ['cluster']] = c
            c += 1

        # assign data samples to clusters again by only using the defining cluster center vectors
        # and find common positions in each cluster
        sample_clustering['nearest_centroid'] = -1

        s = 0
        for sample in enc_2D:
            sample = sample[def_positions]
            name = sample_names[s]

            diff_counts = np.empty(len(def_cluster_snps))
            c = 0
            for center in def_cluster_snps:
                comp = np.equal(sample, center)
                diff_count = len(np.where(comp == False)[0].tolist())
                diff_counts[c] = diff_count
                c += 1

            min_idx = np.argmin(diff_counts, axis=0)

            sample_clustering.loc[sample_clustering['sample_names'] == name, 'nearest_centroid'] = min_idx
            s += 1

        consistent = sample_clustering[
            sample_clustering['cluster'] == sample_clustering['nearest_centroid']].reset_index()
        consist_list = consistent['sample_names'].tolist()

        labels = ["consistent" if x in consist_list else "ambiguous" for x in sample_names]
        umap_df['cluster_assign'] = labels

        fig, ax = plt.subplots(1)
        sns.scatterplot(x='umap_1', y='umap_2', data=umap_df, ax=ax, s=5, hue='cluster_assign')
        sns.scatterplot(x='umap_1', y='umap_2', data=k_Means, ax=ax, s=5, palette='flare',
                        hue=['cluster center'] * len(cluster_centers))
        ax.set_aspect('equal')

        class_list = list(range(0, len(inv_transformed_centers)))
        classes = pd.DataFrame(class_list, columns=['class'])
        classes['kMeans_count'] = classes['class'].map(sample_clustering['cluster'].value_counts())
        classes['check_count'] = classes['class'].map(sample_clustering['nearest_centroid'].value_counts())
        classes = classes.fillna(0)
        # print(classes)

        center_pos = np.where(inv_transformed_points[c] > 0.5)[0].tolist()
        print("Cluster center SNPS: " + str(center_pos))
        print()


def freqGenotype(df, enc_2D, file):
    sample_names = df['label'].tolist()
    cluster_count = df['cluster'].nunique()
    cluster_df = pd.DataFrame(range(cluster_count), columns=['Cluster'])

    # determine SNP frequencies for each cluster
    cluster_snp_freqs = np.full((cluster_count, len(enc_2D[0])), 0.0)

    s = 0
    for sample in enc_2D:
        name = sample_names[s]
        cluster = df[df['label'] == name]['cluster'].values[0]
        cluster_snp_freqs[cluster] = np.add(cluster_snp_freqs[cluster], sample)
        s += 1

    c = 0
    steps = range(5, 11)
    freqs = [[] for _ in steps]
    for cluster in cluster_snp_freqs:
        cluster_size = df[df['cluster'] == c].shape[0]
        cluster_snp_freqs[c] = cluster_snp_freqs[c] / cluster_size

        for step in reversed(steps):
            limit = step / 10
            idx = len(freqs) - 1 - step
            freqs[idx].append(len(np.where(cluster_snp_freqs[c] >= limit)[0]))

        c += 1

    for step in reversed(steps):
        label = str(step / 10)
        idx = len(freqs) - 1 - step
        cluster_df[label] = freqs[idx]

    # create vectors representing each cluster for given threshold
    lim = 0.9
    c = 0
    rep = cluster_snp_freqs.copy()
    for cluster in cluster_snp_freqs:
        rep[c] = np.where(cluster_snp_freqs[c] >= lim, 1, 0)
        c += 1

    cluster_df.to_csv(file, index=False)
    print(cluster_df)
    return rep


def compareToMLST(df, file):
    mlst = MLST.compareMLST("Parr1509_CP004010_SNPSummary.tsv", ["TPANIC_RS00695", "TPANIC_RS02695", "TPANIC_RS03500"],
                            True)
    MLST_df = pd.merge(df, mlst, left_on='label', right_on='Sample', how='inner')

    MLST_grouped = MLST_df.groupby(['umap_cluster'])['Allelic Profile']
    print(MLST_grouped.describe())
    print(df.groupby(['umap_cluster']).size().describe())
    MLST_grouped.describe().to_csv(file)

def randIndex(dfs):
    umap_rand_table = np.empty((len(dfs), len(dfs)))
    tsne_rand_table = np.empty((len(dfs), len(dfs)))
    o = 0
    for df_out in dfs:
        i = 0
        for df_in in dfs:
            rand_clusters = pd.merge(df_out, df_in, on='label', how='inner')
            umap_rand_index = adjusted_rand_score(rand_clusters['umap_cluster_x'].tolist(),
                                                  rand_clusters['umap_cluster_y'].tolist())
            umap_rand_table[o][i] = round(umap_rand_index, 4)
            tsne_rand_index = adjusted_rand_score(rand_clusters['tsne_cluster_x'].tolist(),
                                                  rand_clusters['tsne_cluster_y'].tolist())
            tsne_rand_table[o][i] = round(tsne_rand_index, 4)
            i += 1
        o += 1

    print(umap_rand_table)
    print(tsne_rand_table)


print("No Filter, kMEans, noMLST")
df = compare("Parr1509_CP004010_SNPSummary.tsv", ["LOW", "MODIFIER", "MODERATE", "HIGH"], "binary", "kMeans")
#plotCompare(df, "UMAP vs. tSNE mit k-Means-Clustering", "UMAP", "tSNE")
#compareToMLST(df, 'cluster_allel_profiles.csv')

print("Filter, kMeans, noMLST")
filter_df = compare("Parr1509_CP004010_SNPSummary.tsv", ["MODERATE", "HIGH"], "binary", "kMeans")
#plotCompare(filter_df, "UMAP vs. tSNE mit k-Means-Clustering (MODERATE, HIGH)", "UMAP", "tSNE")
#compareToMLST(filter_df, 'cluster_allel_profiles.csv')

print("No Filter, kMeans, MLST")
MLST_df = compare("Parr1509_CP004010_SNPSummary.tsv", ["LOW", "MODIFIER", "MODERATE", "HIGH"], "binary", "kMeans", ["TPANIC_RS00695", "TPANIC_RS02695", "TPANIC_RS03500"])
#plotCompare(MLST_df, "MLST-UMAP vs. MLST-tSNE mit k-Means-Clustering", "MLST-UMAP", "MLST-tSNE")
#compareToMLST(MLST_df, 'MLST_cluster_allel_profiles.csv')


randIndex([df, filter_df, MLST_df])

####################
"""print("All Data acquired")
##plotCompare(df, "test")


##### MLST ##########
mlst = MLST.compareMLST("Parr1509_CP004010_SNPSummary.tsv", ["TPANIC_RS00695", "TPANIC_RS02695", "TPANIC_RS03500"], True, ["LOW", "MODIFIER", "MEDIUM", "HIGH"])
MLST_df = pd.merge(MLST_df, mlst, left_on='label', right_on='Sample', how='inner')

MLST_df = MLST_df.sort_values(by=['Allelic Profile'])
fig, ax = plt.subplots(1)
fig.subplots_adjust(top=0.8)
sns.scatterplot(x='UMAP_1', y='UMAP_2', data=MLST_df, ax=ax, s=5, hue='Allelic Profile',
                palette=sns.color_palette("hls", MLST_df['Allelic Profile'].nunique()), hue_order=MLST_df['Allelic Profile'].unique().sort())
ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, title='Allel-Profil',
             ncol=legend_cols(MLST_df['Allelic Profile']))
# ax[1].legend([], [], frameon=False)
ax.set_title("MLST UMAP")
plt.tight_layout()
fig.suptitle("MLST UMAP mit Allel-Profilen ", y=0.98)


def MLST_label_UMAP(filter, encoding, cluster_method, loci=[]):
    umap = compare("Parr1509_CP004010_SNPSummary.tsv", filter, encoding, "none")
    MLST_umap = compare("Parr1509_CP004010_SNPSummary.tsv", filter, encoding, cluster_method, loci)
    mlst = MLST.compareMLST("Parr1509_CP004010_SNPSummary.tsv", loci, True, filter)
    print("All Data acquired")

    umap['MLST_cluster'] = MLST_umap['umap_cluster']
    plotCompare(umap, "test")


    fig, ax = plt.subplots(1, 2)
    fig.subplots_adjust(top=0.8)
    sns.scatterplot(x='UMAP_1', y='UMAP_2', data=umap, ax=ax[0], s=5, hue='MLST_cluster',
                    palette=sns.color_palette("hls", umap['MLST_cluster'].nunique()))
    ax[0].set_title("UMAP")
    # ax[0].legend([], [], frameon=False)
    ax[0].legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, title='UMAP Cluster',
                 ncol=legend_cols(umap['MLST_cluster']))

    sns.scatterplot(x='UMAP_1', y='UMAP_2', data=MLST_umap, ax=ax[1], s=5, hue='umap_cluster',
                    palette=sns.color_palette("hls", MLST_umap['umap_cluster'].nunique()))
    ax[1].legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, title='UMAP Cluster',
                 ncol=legend_cols(MLST_umap['umap_cluster']))
    # ax[1].legend([], [], frameon=False)
    ax[1].set_title("UMAP")
    plt.subplots_adjust(wspace=0.5)
    plt.tight_layout()
    fig.suptitle("UMAP with MLST and MLST UMAP labels", y=0.98)"""


# MLST_label_UMAP(["LOW", "MODIFIER", "MEDIUM", "HIGH"], "binary", "kMeans", ["TPANIC_RS00695", "TPANIC_RS02695", "TPANIC_RS03500"])

plt.show()
