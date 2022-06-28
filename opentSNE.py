import os
import sqlite3
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from openTSNE import TSNE


def encodeNuc(nuc):
    enc = ["a", "c", "g", "t", ".", "-"]
    i = enc.index(nuc)
    one_hot = [0] * 6
    one_hot[i] = 1
    return one_hot


class snpTSNE:

    def __init__(self, db_name):
        self.db_name = db_name
        self.symb_numb = 6
        self.default_symb_enc = [0, 0, 0, 0, 1, 0]

        # get tSNE input data from snp database
        # establish connection to sqlite database
        conn = sqlite3.connect(db_name)
        cursor = conn.cursor()

        # get all possible sequence_names aka types
        sql_command = "SELECT DISTINCT sequence_name FROM unambiguous ORDER BY sequence_name;"
        cursor.execute(sql_command)
        content = cursor.fetchall()
        self.sequence_names = list()
        for seq in content:
            self.sequence_names.append(seq[0])
        self.seq_count = len(self.sequence_names)

        # get all reference_GenWidePos for SNPs
        sql_command = "SELECT DISTINCT reference_GenWidePos FROM unambiguous ORDER BY reference_GenWidePos;"
        cursor.execute(sql_command)
        content = cursor.fetchall()
        self.positions = list()
        for pos in content:
            self.positions.append(pos[0])
        self.pos_count = len(self.positions)

        # create training vector of SNPs, with field for each SNP site (reference_GenWidePos) using one-hot encoding
        reference_3D = np.full((1, self.pos_count, self.symb_numb), self.default_symb_enc)
        train_3D = np.full((self.seq_count, self.pos_count, self.symb_numb), self.default_symb_enc)

        for sequence in self.sequence_names:
            # query SNP pattern and position from database
            sql_command = "SELECT SNP_pattern, reference_GenWidePos FROM unambiguous WHERE unambiguous.sequence_name='" + str(sequence) + "' ORDER BY reference_GenWidePos;"
            cursor.execute(sql_command)
            snp_query = cursor.fetchall()

            for entry in snp_query:
                pos = self.positions.index(entry[1])
                seq = self.sequence_names.index(sequence)
                enc_entry = encodeNuc(entry[0][0])
                train_3D[seq][pos] = enc_entry

                reference_3D[0][pos] = encodeNuc(entry[0][1])
        conn.close()

        train_3D=np.append(train_3D, reference_3D, axis=0)

        # reshape training vector from 3D to 2D
        sample, position, one_hot_enc = train_3D.shape
        train_2D = train_3D.reshape((sample,position*one_hot_enc))

        self.tsne = TSNE(
            perplexity=30,
            metric="euclidean",
            n_jobs=8,
            random_state=42,
            verbose=True,
        )
        self.embedding_train = self.tsne.fit(train_2D)
        self.labels = self.sequence_names.copy()
        self.labels.append("Reference(NC_021490)")


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

            print(tsne_result_df["Strain"])
            plt.show()




    def testDataTSNE(self, tsvFile):
        if os.path.isfile(tsvFile):
            query = pd.read_csv(tsvFile, sep='\t', header=1)
            # create vector of SNPs, with field for each SNP site (reference_GenWidePos)
            sample_names = list(query.columns)
            sample_names.remove('Position')
            # sample_names.remove('Reference')

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

            #embedding_eval = self.embedding_train.transform(eval_2D)
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
            plt.show()



test = snpTSNE("snps.db")
#test.embedQuery("PT_SIF0908.variants.tsv")
#test.multembedTSNE("variantContentTable.tsv")
#test.testDataTSNE("variantContentTable.tsv")

def adaptTSNE():
    t = 1