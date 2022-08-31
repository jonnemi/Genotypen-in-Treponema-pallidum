import os
import sqlite3
import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import seaborn as sns

from openTSNE import TSNE
import dataProcess

def adaTSNE(queryTSV, db_name):
    dataset = dataProcess.getADAdataset(queryTSV, db_name, "one-hot")
    strains = dataProcess.assignStrains(dataset['test_seq_names'])

    # first tSNE using only query data
    tsne1 = TSNE(
        perplexity=30,
        initialization="pca",
        metric="manhattan",
        n_jobs=8,
        random_state=3,
        #verbose=True,
    )
    reducer1 = tsne1.fit(dataset["test"])
    only_query = reducer1.transform(dataset["test"])

    only_query_df = pd.DataFrame(
        {'tSNE_1': only_query[:, 0], 'tSNE_2': only_query[:, 1],
         'label': dataset['test_seq_names']})
    only_query_df["strain"] = strains

    # second tSNE trained on reference genomes with embedded query data
    tsne2 = TSNE(
        perplexity=30,
        initialization="pca",
        metric="manhattan",
        n_jobs=8,
        random_state=3,
        #verbose=True,
    )
    reducer2 = tsne2.fit(dataset["train"])
    train = reducer2.transform(dataset["train"])
    test = reducer2.transform(dataset["test"])
    train_df = pd.DataFrame(
        {'tSNE_1': train[:, 0], 'tSNE_2': train[:, 1],
         'label': dataset['train_seq_names']})
    train_df["strain"] = "Referenz"

    test_df = pd.DataFrame(
        {'tSNE_1': test[:, 0], 'tSNE_2': test[:, 1],
         'label': dataset['test_seq_names']})
    test_df["strain"] = strains

    with_ref_df = train_df.append(test_df)

    # plot both tSNEs side by side
    fig, ax = plt.subplots(2, 1)
    fig.subplots_adjust(top=0.8)
    sns.scatterplot(x='tSNE_1', y='tSNE_2', hue='strain', data=only_query_df, ax=ax[0], s=20, style='strain',
                    markers=["D", "v"], palette=["C0", "C1"])
    ax[0].set_title("Nur Query-Daten")
    ax[0].legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, title='Stamm')

    sns.scatterplot(x='tSNE_1', y='tSNE_2', hue='strain', data=with_ref_df, ax=ax[1], s=20, style='strain',
                    markers=["o", "D", "v"], palette=["k", "C0", "C1"])
    ax[1].legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, title='Stamm')
    # ax[1].legend([], [], frameon=False)
    ax[1].set_title("Referenz-tSNE mit Query-Daten")
    # ax[1].set_aspect('equal')
    plt.subplots_adjust(wspace=0.5)
    plt.tight_layout()
    fig.suptitle("tSNE Vergleich", y=0.98)
    #plt.show()

#adaTSNE("66-Proben-Datensatz.tsv", "snps.db")

def encodeNuc(nuc):
    enc = ["a", "c", "g", "t", ".", "-"]
    i = enc.index(nuc)
    one_hot = [0] * 6
    one_hot[i] = 1
    return one_hot


class snpTSNE:

    def __init__(self, db_name):
        self.db_name = db_name
        self.default_symb_enc = [0, 0, 0, 0, 1]

        self.SNPdf = dataProcess.getSNPdf(db_name)
        self.ref_positions = self.SNPdf["Position"].unique().tolist()
        self.ref_positions.sort()
        self.sequence_names = self.SNPdf['Sequence_name'].unique().tolist()
        self.sequence_names.sort()
        self.sequence_names.append("Reference(NC_021490)")
        train_2D = dataProcess.encDF(self.SNPdf, self.default_symb_enc, self.ref_positions)

        """self.tsne = TSNE(
            perplexity=30,
            metric="euclidean",
            n_jobs=8,
            random_state=42,
            verbose=True,
        )"""
        self.tsne = TSNE(
            perplexity=30,
            initialization="random",
            metric="euclidean",
            n_jobs=8,
            random_state=3,
            verbose=True
        )
        self.embedding_train = self.tsne.fit(train_2D)

    def embedQuery(self, tsvFile):
        if os.path.isfile(tsvFile):
            query = pd.read_csv(tsvFile, sep='\t', header=None, skiprows=1).to_numpy()

            # create vector of query SNPs, with field for each SNP site (reference_GenWidePos) using one-hot encoding
            test_3D = np.full((1, self.pos_count, self.symb_numb), self.default_symb_enc)
            for entry in query:
                if len(entry[1]) == 1 and entry[0] in self.positions:
                    pos = self.positions.index(entry[0])
                    enc_entry = encodeNuc(entry[1][0].lower())
                    test_3D[0][pos] = enc_entry

            # reshape training vector from 3D to 2D
            sample, position, one_hot_enc = test_3D.shape
            test_2D = test_3D.reshape((sample, position * one_hot_enc))

            embedding_test = self.embedding_train.transform(test_2D)

            # plot query embedded in train tSNE
            train_tsne_df = pd.DataFrame(
                {'tsne_1': self.embedding_train[:, 0], 'tsne_2': self.embedding_train[:, 1],
                 'label': self.labels})
            test_tsne_df = pd.DataFrame(
                {'tsne_1': embedding_test[:, 0], 'tsne_2': embedding_test[:, 1],
                 'label': ["Query"]})
            fig, ax = plt.subplots(1)
            sns.scatterplot(x='tsne_1', y='tsne_2', hue='label', data=train_tsne_df, ax=ax, s=5, alpha=0.5)
            sns.scatterplot(x='tsne_1', y='tsne_2', hue='label', data=test_tsne_df, ax=ax, s=5, alpha=1)
            ax.set_aspect('equal')
            ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0, ncol=3)
            fig.tight_layout()
            plt.show()

            # return list of Euclidean Distances to new embedded test point
            eukl_dists = list()
            for point in self.embedding_train:
                eukl_dists.append(np.linalg.norm(np.array(point) - np.array(embedding_test[0])))

            min_index = np.argmin(eukl_dists)
            nearest = self.sequence_names[min_index]
            print(nearest)
            return (eukl_dists)

    def multembedTSNE(self, tsvFile):
        if os.path.isfile(tsvFile):
            query = pd.read_csv(tsvFile, sep='\t', header=1)
            # create vector of SNPs, with field for each SNP site (reference_GenWidePos)
            sample_names = list(query.columns)
            sample_names.remove('Position')
            sample_names.remove('Reference')

            # remove all rows where position is not in reference database
            query = query[query['Position'].isin(map(str, self.positions))]
            query.reset_index()

            # create vector query SNPs, with field for each SNP site (reference_GenWidePos) using one-hot encoding
            eval_3D = np.full((len(sample_names), len(self.positions), self.symb_numb), self.default_symb_enc)

            s = 0
            for sample in sample_names:
                p = 0
                for entry in query[sample]:
                    pos = self.positions.index(int(query['Position'].iloc[p]))
                    eval_3D[s][pos] = encodeNuc(entry.lower())
                    p += 1
                s += 1

            # reshape training vector from 3D to 2D
            sample, position, one_hot_enc = eval_3D.shape
            eval_2D = eval_3D.reshape((sample, position * one_hot_enc))

            embedding_eval = self.embedding_train.transform(eval_2D)

            group1 = ['BosniaA']
            group2 = ['Fribourg', 'CDC2', 'GHA1', 'Gauthier', 'IND1', 'SAM1', 'SamoaD']
            group3 = ['Seattle81', 'SEA86', 'NE20', 'BAL3', 'BAL73', 'Chicago', 'NIC1', 'NIC2', 'Dallas']

            strain_labels = ['TEN', 'TPE', 'Nicols (TPA)', 'SS14 (TPA)']

            strains = []
            for sample in sample_names:
                if sample in group1:
                    strain = strain_labels[0]
                elif sample in group2:
                    strain = strain_labels[1]
                elif sample in group3:
                    strain = strain_labels[2]
                else:
                    strain = strain_labels[3]
                strains.append(strain)

            train_tsne_df = pd.DataFrame(
                {'tsne_1': self.embedding_train[:, 0], 'tsne_2': self.embedding_train[:, 1],
                 'label': self.labels})
            tsne_result_df = pd.DataFrame(
                {'tsne_1': embedding_eval[:, 0], 'tsne_2': embedding_eval[:, 1], 'Sample': sample_names})
            tsne_result_df["Strain"] = strains
            fig, ax = plt.subplots(1)
            sns.scatterplot(x='tsne_1', y='tsne_2', data=train_tsne_df, ax=ax, s=5, alpha=0.5, color="black")
            sns.scatterplot(x='tsne_1', y='tsne_2', hue='Sample', data=tsne_result_df, ax=ax, s=5, style='Strain',
                            markers=["o", "v", "D", "X"])
            ax.set_aspect('equal')
            ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0, ncol=3)
            fig.tight_layout()
            plt.show()



    def adaptTSNE(self, tsvFile):
        if os.path.isfile(tsvFile):
            query = pd.read_csv(tsvFile, sep='\t', header=1)

        # manually drop TPE an TEN samples
        query = query.drop(columns=['BosniaA', 'Fribourg', 'CDC2', 'GHA1', 'Gauthier', 'IND1', 'SAM1', 'SamoaD'])

        # create vector of SNPs, with field for each SNP site (reference_GenWidePos)
        sample_names = list(query.columns)
        sample_names.remove('Position')
        sample_names.remove('Reference')

        # find common SNP positions in reference database and query
        query_pos = query.Position.tolist()
        common_pos = list(set(map(str, self.ref_positions)).intersection(query_pos))
        common_pos.sort(key=lambda pos: int(pos))

        # remove all rows from dataframes that don't have common reference genome position for SNP
        query = query[query['Position'].isin(common_pos)]
        query['Position'] = pd.to_numeric(query['Position'])
        query.reset_index()
        common_pos = list(map(int, common_pos))
        SNPdf = self.SNPdf[self.SNPdf['Position'].isin(common_pos)]
        SNPdf.reset_index()

        # create training vector of SNPs, with field for each SNP site (reference_GenWidePos) using one-hot encoding
        train_2D = dataProcess.encDF(self.SNPdf, self.default_symb_enc, common_pos)
        embedding_ADAtrain = self.tsne.fit(train_2D)

        # create vector of query SNPs, with field for each SNP site (reference_GenWidePos) using one-hot encoding
        eval_2D = dataProcess.encQuery(query, self.default_symb_enc, common_pos)
        embedding_eval = embedding_ADAtrain.transform(eval_2D)

        # assign strains to query data
        strains = dataProcess.assignStrains(sample_names)

        # plot tSNE
        train_tsne_df = pd.DataFrame(
            {'tsne_1': embedding_ADAtrain[:, 0], 'tsne_2': embedding_ADAtrain[:, 1],
             'label': self.sequence_names, 'strain': "unknown (Reference db)"})
        tsne_result_df = pd.DataFrame(
            {'tsne_1': embedding_eval[:, 0], 'tsne_2': embedding_eval[:, 1], 'label': sample_names, 'strain': strains})
        fig, ax = plt.subplots(1)
        sns.scatterplot(x='tsne_1', y='tsne_2', data=train_tsne_df, ax=ax, s=5, alpha=0.5, color="black")
        sns.scatterplot(x='tsne_1', y='tsne_2', hue='label', data=tsne_result_df, ax=ax, s=5, style='strain',
                        markers=["o", "v"])
        ax.set_aspect('equal')
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0, ncol=3)
        fig.tight_layout()

        # save tSNE dataframes to file
        tSNE_df = pd.concat([train_tsne_df, tsne_result_df], ignore_index=True)
        tSNE_df.to_csv("tSNE.tsv", sep='\t')

        plt.show()


#
# test.embedQuery("PT_SIF0908.variants.tsv")
# test.multembedTSNE("66-Proben-Datensatz.tsv")
# test.testDataTSNE("66-Proben-Datensatz.tsv")
#test.adaptTSNE("66-Proben-Datensatz.tsv")

def MLSTtSNE(tsvFile):
    dataset = dataProcess.getLociDataset(tsvFile, "snps.db", [0, 0, 0, 0, 1],
                                         ("TPANIC_RS00695", "TPANIC_RS02695", "TPANIC_RS03500"))
    print(dataset["train"])

    embedding = TSNE(
        perplexity=23,
        initialization="pca",
        metric="cosine",
        n_jobs=8,
        random_state=3,

    ).fit(dataset["test"])

    # assign strains to query data
    strains = dataProcess.assignStrains(dataset['test_sample_names'])

    # plot tSNE
    """train_tsne_df = pd.DataFrame(
        {'tsne_1': fit[:, 0], 'tsne_2': fit[:, 1],
         'label': dataset['train_seq_names'], 'strain': "unknown (Reference db)"})"""
    tsne_result_df = pd.DataFrame(
        {'tsne_1': embedding[:, 0], 'tsne_2': embedding[:, 1], 'label': dataset['test_sample_names'],
         'strain': strains})
    fig, ax = plt.subplots(1)
    #sns.scatterplot(x='tsne_1', y='tsne_2', data=train_tsne_df, ax=ax, s=5, alpha=0.5, color="black")
    sns.scatterplot(x='tsne_1', y='tsne_2', hue='strain', data=tsne_result_df, ax=ax, s=5, style='strain',
                    markers=["o", "v"])
    ax.set_aspect('equal')
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0, ncol=3)
    fig.tight_layout()

    # save MLST tSNE dataframes to file
    # tSNE_df = pd.concat([train_tsne_df, tsne_result_df], ignore_index=True)
    # tSNE_df.to_csv("MLSTtSNE.tsv", sep='\t')

    plt.show()

#MLSTtSNE("66-Proben-Datensatz.tsv")


def testDataTSNE(tsvFile, default_enc):
    if os.path.isfile(tsvFile):
        query = pd.read_csv(tsvFile, sep='\t', header=1)
        # create vector of SNPs, with field for each SNP site (reference_GenWidePos)
        sample_names = list(query.columns)
        sample_names.remove('Position')
        # sample_names.remove('Reference')

        # remove all rows with insertions
        query = query[~query['Position'].str.contains("\+")]
        query.reset_index()


        # create vector query SNPs, with field for each SNP site (reference_GenWidePos) using one-hot encoding
        eval_3D = np.empty(len(sample_names), int(len(query['Position'])), len(default_enc))

        s = 0
        for sample in sample_names:
            p = 0
            for entry in query[sample]:
                pos = int(query['Position'].iloc[p])
                eval_3D[s][pos] = encodeNuc(entry.lower())
                p += 1
            s += 1

        # reshape training vector from 3D to 2D
        sample, position, one_hot_enc = eval_3D.shape
        eval_2D = eval_3D.reshape((sample, position * one_hot_enc))
        print(eval_2D)

        """# embedding_eval = self.embedding_train.transform(eval_2D)
        embedding_eval = self.tsne.fit(eval_2D)

        # plot training tSNE

        # assign color palette to each sample according to it's expected group
        color1 = "pink"
        color2 = "blue"
        color3 = "green"
        color4 = "orange"

        group1 = ['BosniaA']
        group2 = ['Fribourg', 'CDC2', 'GHA1', 'Gauthier', 'IND1', 'SAM1', 'SamoaD']
        group3 = ['Seattle81', 'SEA86', 'NE20', 'BAL3', 'BAL73', 'Chicago', 'NIC1', 'NIC2', 'Dallas']

        strain_labels = ['TEN', 'TPE', 'Nicols (TPA)', 'SS14 (TPA)']

        palette = dict()
        strains = []
        for sample in sample_names:
            if sample == 'Reference':
                color = 'red'
                strain = strain_labels[2]
            elif sample in group1:
                color = color1
                strain = strain_labels[0]
            elif sample in group2:
                color = color2
                strain = strain_labels[1]
            elif sample in group3:
                color = color3
                strain = strain_labels[2]
            else:
                color = color4
                strain = strain_labels[3]
            palette[sample] = color
            strains.append(strain)

        tsne_result_df = pd.DataFrame(
            {'tsne_1': embedding_eval[:, 0], 'tsne_2': embedding_eval[:, 1], 'Sample': sample_names})
        tsne_result_df["Strain"] = strains
        fig, ax = plt.subplots(1)
        sns.scatterplot(x='tsne_1', y='tsne_2', hue='Sample', data=tsne_result_df, ax=ax, s=5, style='Strain',
                        markers=["o", "v", "D", "X"])
        ax.set_aspect('equal')
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0, ncol=3)
        fig.tight_layout()
        plt.show()"""

#testDataTSNE("66-Proben-Datensatz.tsv", [0, 0, 0, 0, 1])