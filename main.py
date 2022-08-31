import database
import MLST
import UMAP
import opentSNE
import tSNEvsUMAP
from Bio import SeqIO
import sys, os
import matplotlib.pyplot as plt
from datetime import datetime
startTime = datetime.now()


# Disable print
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

# Restore print
def enablePrint():
    sys.stdout = sys.__stdout__


if __name__ == '__main__':
    seperator_len = 80
    sep = "=" * seperator_len
    sub_sep = "-" * seperator_len
    MLST_loci = ["TPANIC_RS00695", "TPANIC_RS02695", "TPANIC_RS03500"]

    # create reference genome database snps.db
    database.loadDB()

    # reference genome
    genome_record = SeqIO.read("NC_021490.2.gb", "genbank")

    # MLST
    print(sep)
    MLST.lociSNVdensity(MLST.getRefLoci(genome_record))
    print("MLST des 66-Proben-Datensatzes:")
    MLST.compareMLST("66-Proben-Datensatz.tsv", MLST_loci, False)
    print(sub_sep)
    print("MLST des 1508-Proben-Datensatzes ohne Filtern nach SNP-Auswirkungen:")
    MLST.compareMLST("1508-Proben-Datensatz.tsv", MLST_loci, True)
    print(sep)
    print("MLST des 1508-Proben-Datensatzes mit Filtern nach MODERATE und HIGH SNP-Auswirkungen:")
    MLST.compareMLST("1508-Proben-Datensatz.tsv", MLST_loci, True, ["MODERATE", "HIGH"])
    print(sep)

    # UMAP 66-Proben-Datensatz
    print("UMAP des 66-Proben-Datensatzes:")
    UMAP.adaUmap("66-Proben-Datensatz.tsv", "snps.db")
    print(sep)

    # tSNE 66-Proben-Datensatz
    print("tSNE des 66-Proben-Datensatzes:")
    opentSNE.adaTSNE("66-Proben-Datensatz.tsv", "snps.db")
    print(sep)

    # tSNE vs. UMAP für 1508-Proben-Datensatz
    print("tSNE vs. UMAP für den 1508-Proben-Datensatz (ohne SNP-Filtern, ohne MLST)")
    print(sub_sep)
    df = tSNEvsUMAP.compare("1508-Proben-Datensatz.tsv", ["LOW", "MODIFIER", "MODERATE", "HIGH"], "binary", "kMeans")
    tSNEvsUMAP.plotCompare(df, "UMAP vs. tSNE mit k-Means-Clustering", "UMAP", "tSNE")
    blockPrint()
    cluster_profiles = tSNEvsUMAP.compareToMLST(df, 'cluster_allel_profiles.csv')
    enablePrint()
    print("Allel-Profile der UMAP-Cluster:")
    print(cluster_profiles)
    print(sep)

    print("tSNE vs. UMAP für den 1508-Proben-Datensatz (mit SNP-Filtern nach MODERATE + HIGH, ohne MLST):")
    print(sub_sep)
    filter_df = tSNEvsUMAP.compare("1508-Proben-Datensatz.tsv", ["MODERATE", "HIGH"], "binary", "kMeans")
    tSNEvsUMAP.plotCompare(filter_df, "UMAP vs. tSNE mit k-Means-Clustering (MODERATE, HIGH)", "UMAP", "tSNE")
    blockPrint()
    tSNEvsUMAP.compareToMLST(filter_df, 'cluster_allel_profiles.csv')
    enablePrint()
    print("Allel-Profile der UMAP-Cluster:")
    print(cluster_profiles)
    print(sep)

    print("tSNE vs. UMAP für den 1508-Proben-Datensatz (ohne SNP-Filtern, mit MLST):")
    print(sub_sep)
    MLST_df = tSNEvsUMAP.compare("1508-Proben-Datensatz.tsv", ["LOW", "MODIFIER", "MODERATE", "HIGH"], "binary", "kMeans", MLST_loci)
    tSNEvsUMAP.plotCompare(MLST_df, "MLST-UMAP vs. MLST-tSNE mit k-Means-Clustering", "MLST-UMAP", "MLST-tSNE")
    blockPrint()
    tSNEvsUMAP.compareToMLST(MLST_df, 'MLST_cluster_allel_profiles.csv')
    enablePrint()
    print("Allel-Profile der UMAP-Cluster:")
    print(cluster_profiles)
    print(sep)

    print("Adjusted-Rand-Index der berechneten tSNEs und UMAPs:")
    tSNEvsUMAP.randIndex([df, filter_df, MLST_df])

print()
print(sep)
print("Laufzeit: " + str(datetime.now() - startTime))
plt.show()