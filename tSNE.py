import numpy as np
import scipy.cluster as sc_clstr
import sqlite3
from numpy import array
from numpy import argmax
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder

from sklearn.manifold import TSNE
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def encodeNuc(nuc):
    enc = [".", "a", "c", "g", "t"]
    i = enc.index(nuc)
    one_hot = [0] * 5
    one_hot[i] = 1
    return i

#def generate_sample_vector():
    #return np.array([float(n) for n in np.random.randint(1, 4, size=vector_size)])
# sample_data = np.array([generate_sample_vector() for i in range(samples)])

# get tSNE input data from snp database
# establish connection to sqlite database
conn = sqlite3.connect("snps.db")
cursor = conn.cursor()

# get all possible sequence_names aka types
sql_command = "SELECT DISTINCT sequence_name FROM unambiguous ORDER BY sequence_name;"
cursor.execute(sql_command)
content = cursor.fetchall()
sequence_names = list()
for seq in content:
    sequence_names.append(seq[0])

# get all reference_GenWidePos for SNPs
sql_command = "SELECT DISTINCT reference_GenWidePos FROM unambiguous ORDER BY reference_GenWidePos;"
cursor.execute(sql_command)
content = cursor.fetchall()
positions = list()
for pos in content:
    positions.append(pos[0])


# create vector of SNPs, with field for each SNP site (reference_GenWidePos)
"""vec = np.full((len(sequence_names), len(positions), 5), [1, 0, 0, 0, 0])
for sequence in sequence_names:
    # query SNP pattern and position from database
    sql_command = "SELECT SNP_pattern, reference_GenWidePos FROM unambiguous WHERE unambiguous.sequence_name='" + str(sequence) + "' ORDER BY reference_GenWidePos;"
    cursor.execute(sql_command)
    snp_query = cursor.fetchall()

    for entry in snp_query:
        pos = positions.index(entry[1])
        seq = sequence_names.index(sequence)
        enc_entry = encodeNuc(entry[0][0])
        vec[seq][pos] = enc_entry"""

### 2nd try with dim 2
vec = np.full((len(sequence_names), len(positions)), 0)
for sequence in sequence_names:
    # query SNP pattern and position from database
    sql_command = "SELECT SNP_pattern, reference_GenWidePos FROM unambiguous WHERE unambiguous.sequence_name='" + str(sequence) + "' ORDER BY reference_GenWidePos;"
    cursor.execute(sql_command)
    snp_query = cursor.fetchall()

    for entry in snp_query:
        pos = positions.index(entry[1])
        seq = sequence_names.index(sequence)
        enc_entry = encodeNuc(entry[0][0])
        vec[seq][pos] = enc_entry
conn.close()

tsne = TSNE(2)
tsne_result = tsne.fit_transform(vec)
print(tsne_result.shape)
tsne_result_df = pd.DataFrame({'tsne_1': tsne_result[:, 0], 'tsne_2': tsne_result[:, 1], 'label': sequence_names})
fig, ax = plt.subplots(1)
sns.scatterplot(x='tsne_1', y='tsne_2', hue='label', data=tsne_result_df, ax=ax, s=5)
# lim = (tsne_result.min()-5, tsne_result.max()+5)
# ax.set_xlim(lim)
# ax.set_ylim(lim)
ax.set_aspect('equal')
ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)
plt.show()
