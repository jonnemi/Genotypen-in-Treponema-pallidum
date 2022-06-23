import math
import os
import re
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
    return i


class snpTSNE:

    def __init__(self, db_name):
        self.db_name = db_name

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

        # create vector of SNPs, with field for each SNP site (reference_GenWidePos)
        train = np.full((self.seq_count, self.pos_count), 0)
        for sequence in self.sequence_names:
            # query SNP pattern and position from database
            sql_command = "SELECT SNP_pattern, reference_GenWidePos FROM unambiguous WHERE unambiguous.sequence_name='" + str(
                sequence) + "' ORDER BY reference_GenWidePos;"
            cursor.execute(sql_command)
            snp_query = cursor.fetchall()

            for entry in snp_query:
                pos = self.positions.index(entry[1])
                seq = self.sequence_names.index(sequence)
                enc_entry = encodeNuc(entry[0][0])
                train[seq][pos] = enc_entry
        conn.close()

        self.tsne = TSNE(
            perplexity=30,
            metric="euclidean",
            n_jobs=8,
            random_state=42,
            verbose=True,
        )
        self.embedding_train = self.tsne.fit(train)

        """tsne_result_df = pd.DataFrame(
            {'tsne_1': self.embedding_train[:, 0], 'tsne_2': self.embedding_train[:, 1], 'label': self.sequence_names})
        fig, ax = plt.subplots(1)
        sns.scatterplot(x='tsne_1', y='tsne_2', hue='label', data=tsne_result_df, ax=ax, s=5)
        ax.set_aspect('equal')
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)
        plt.show()"""

    def embedQuery(self, tsvFile):
        if os.path.isfile(tsvFile):
            query = pd.read_csv(tsvFile, sep='\t', header=None, skiprows=1).to_numpy()
            # create vector of SNPs, with field for each SNP site (reference_GenWidePos)
            test = np.full((1, self.pos_count), 0)
            for entry in query:
                if len(entry[1]) == 1 and entry[0] in self.positions:
                    pos = self.positions.index(entry[0])
                    enc_entry = encodeNuc(entry[1][0].lower())
                    test[0][pos] = enc_entry
            embedding_test = self.embedding_train.transform(test)

            train_tsne_df = pd.DataFrame(
                {'tsne_1': self.embedding_train[:, 0], 'tsne_2': self.embedding_train[:, 1],
                 'label': self.sequence_names})
            test_tsne_df = pd.DataFrame(
                {'tsne_1': embedding_test[:, 0], 'tsne_2': embedding_test[:, 1],
                 'label': ["query"]})
            fig, ax = plt.subplots(1)
            sns.scatterplot(x='tsne_1', y='tsne_2', hue='label', data=train_tsne_df, ax=ax, s=5, alpha=0.5)
            sns.scatterplot(x='tsne_1', y='tsne_2', hue='label', data=test_tsne_df, ax=ax, s=5, alpha=1)
            ax.set_aspect('equal')
            ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)
            plt.show()

            # return list of Euclidean Distances to new embedded test point
            eukl_dists = list()
            for point in self.embedding_train:
                eukl_dists.append(np.linalg.norm(np.array(point) - np.array(embedding_test[0])))

            return(eukl_dists)

    def evaluateTSNE(self, tsvFile):
        if os.path.isfile(tsvFile):
            query = pd.read_csv(tsvFile, sep='\t', header=1)
            # create vector of SNPs, with field for each SNP site (reference_GenWidePos)
            sample_names = list(query.columns)
            sample_names.remove('Position')
            sample_names.remove('Reference')

            # remove all rows where position is not in reference database
            query = query[query['Position'].isin(map(str, self.positions))]
            query.reset_index()


            for sample in query[sample_names]:
                vec = np.full((1, self.pos_count), 0)
                for i in range(0, len(query[sample]) - 1):
                    snp = query[sample][i]
                    pos = query['Position'][i]
                    #j = self.positions.index(pos)
                    #vec[0][j] = encodeNuc(snp)
                print(vec)
                """[0] in self.positions:
                    pos = self.positions.index(entry[0])
                    enc_entry = encodeNuc(entry[1][0].lower())
                    test[0][pos] = enc_entry
        min_index = np.argmin(dists)
        nearest = self.sequence_names[min_index]
        print(nearest)

        print(type(pd.read_csv(tsvFile).columns[0]))

        real_type = "NZ_CP016049"
        print(dists[self.sequence_names.index(real_type)])"""

test = snpTSNE("snps.db")
# test.embedQuery("PT_SIF0908.variants.tsv")
test.evaluateTSNE("variantContentTable.tsv")
