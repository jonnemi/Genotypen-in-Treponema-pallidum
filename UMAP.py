import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import dataProcess
import umap
import os
import kMeans
import dbscan
from collections import Counter
import ast
from scipy.spatial.distance import cdist


def adaUmap(queryTSV, db_name):
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

#adaUmap("variantContentTable.tsv", "snps.db")


def queryUmap(tsvFile, filter, default_enc, loci=list()):
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
    seen = set()
    dupes = [x for x in positions if x in seen or seen.add(x)]
    #print("Duplicate SNP positions: " + str(dupes))

    # check if MLST is wanted
    if loci:
        # load loci dataframe
        loci_df = pd.read_csv("loci.tsv", sep='\t', converters={"var sites": ast.literal_eval})
        # create list of MLST loci SNP positions
        MLST_loci = loci_df[loci_df['locus_tag'].isin(loci)]
        MLST_positions = []
        MLST_loci["var sites"].apply(lambda x: MLST_positions.extend(x))

        rRNA_R8_v1 = 235204
        rRNA_R9_v1 = 235205
        rRNA_R8_v2 = 283649
        rRNA_R9_v2 = 283650

        rRNA_pos = [rRNA_R8_v1, rRNA_R9_v1, rRNA_R8_v2, rRNA_R9_v2]

        MLST_positions.extend(rRNA_pos)
        MLST_positions.sort()

        query = query[query["Position"].isin(MLST_positions)].reset_index()


    #dataProcess.SNPdf_toMega(query, "all_moderate_high.meg")

    # create encoded SNP vector, with field for each SNP site using the given encoding
    process = dataProcess.newEncQuery(query, default_enc, positions)
    enc_2D = process["data"]
    sample_names = process["sample_names"]

    print("Dataset acquired, starting UMAP...")

    # apply UMAP to data
    reducer = umap.UMAP(metric='manhattan',
                        n_neighbors = 30,
                        min_dist = 0.0,
                        n_components = 2,
                        random_state = 42,
                        verbose=True
                        )
    transform = reducer.fit_transform(enc_2D)

    # process umap dimension reduction with kMeans
    umap_df = pd.DataFrame(
        {'umap_1': transform[:, 0], 'umap_2': transform[:, 1],
         'label': sample_names})


    # define clusters in UMAP projection using k-means
    k_Means = kMeans.UMAP_KMeans(umap_df)

    # inverse transform kMeans cluster centers into high dimensional space
    cluster_centers = k_Means[['umap_1', 'umap_2']].to_numpy()
    inv_transformed_points = reducer.inverse_transform(cluster_centers)

    # tune threshold that bins continuous features of high dimensional cluster centers into binary
    tune_threshold = [0.001, 0.001, 0.01]
    tune_threshold.extend(np.linspace(0,1,11))

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


        consistent = sample_clustering[sample_clustering['cluster'] == sample_clustering['nearest_centroid']].reset_index()
        consist_list = consistent['sample_names'].tolist()

        labels = ["consistent" if x in consist_list else "ambiguous" for x in sample_names]
        umap_df['cluster_assign'] = labels

        #fig, ax = plt.subplots(1)
        #sns.scatterplot(x='umap_1', y='umap_2', data=umap_df, ax=ax, s=5, hue='cluster_assign')
        #sns.scatterplot(x='umap_1', y='umap_2', data=k_Means, ax=ax, s=5, palette='flare', hue=['cluster center']*len(cluster_centers))
        #ax.set_aspect('equal')


        class_list = list(range(0, len(inv_transformed_centers)))
        classes = pd.DataFrame(class_list, columns=['class'])
        classes['kMeans_count'] = classes['class'].map(sample_clustering['cluster'].value_counts())
        classes['check_count'] = classes['class'].map(sample_clustering['nearest_centroid'].value_counts())
        classes = classes.fillna(0)
        #print(classes)

    # determine SNP frequencies for each cluster
    cluster_snp_freqs = np.full((len(cluster_centers), len(positions)), 0.0)

    s = 0
    for sample in enc_2D:
        name = sample_names[s]
        cluster = sample_clustering[sample_clustering['sample_names'] == name]['cluster'].values[0]
        cluster_snp_freqs[cluster] = np.add(cluster_snp_freqs[cluster], sample)
        s += 1

    c = 0

    print("compare defining snps")
    print(inv_transformed_points)

    for cluster in cluster_snp_freqs:
        print("Cluster " + str(c))
        cluster_size = k_Means[k_Means['K-means_cluster'] == c]['size'].values[0]
        cluster_snp_freqs[c] = cluster_snp_freqs[c] / cluster_size

        all = np.where(cluster_snp_freqs[c] == 1)[0].tolist()
        print("SNPs occuring in all sequences in this cluster: " + str(all))

        most = np.where(cluster_snp_freqs[c] >= 0.5)[0].tolist()
        print("SNPs occuring in >= 50% of sequences in this cluster: " + str(most))


        center_pos = np.where(inv_transformed_points[c] > 0.5)[0].tolist()
        print("Cluster center SNPS: " + str(center_pos))
        print()
        c += 1

    plt.show()


queryUmap("Parr1509_CP004010_SNPSummary.tsv", ["MODERATE", "HIGH"], "binary", ["TPANIC_RS00695", "TPANIC_RS02695", "TPANIC_RS03500"])

