import database
import dataProcess
import MLST
import UMAP
import opentSNE
import tSNEvsUMAP
from Bio import SeqIO
import matplotlib.pyplot as plt


if __name__ == '__main__':
    sep = "=" * 30
    sub_sep = "-" * 30
    # reference genome
    genome_record = SeqIO.read("NC_021490.2.gb", "genbank")

    # MLST
    print(sep)
    MLST.lociSNVdensity(MLST.getRefLoci(genome_record))
    print("MLST des 66-Proben-Datensatzes:")
    MLST.compareMLST("variantContentTable.tsv", ["TPANIC_RS00695", "TPANIC_RS02695", "TPANIC_RS03500"], False)
    print(sub_sep)
    print("MLST des 1508-Proben-Datensatzes:")
    MLST.compareMLST("Parr1509_CP004010_SNPSummary.tsv", ["TPANIC_RS00695", "TPANIC_RS02695", "TPANIC_RS03500"], True, ["MODERATE", "HIGH"])
    print(sep)

    # UMAP 66-Proben-Datensatz
    print("UMAP des 66-Proben-Datensatzes:")
    UMAP.adaUmap("variantContentTable.tsv", "snps.db")
    print(sep)

    # tSNE 66-Proben-Datensatz
    print("tSNE des 66-Proben-Datensatzes:")
    opentSNE.adaTSNE("variantContentTable.tsv", "snps.db")
    print(sep)

    # tSNE vs. UMAP für 1508-Proben-Datensatz
    print("tSNE vs. UMAP für den 1508-Proben-Datensatz:")
    print(sep)
    print("Ohne SNP-Filtern, ohne MLST:")
    df = tSNEvsUMAP.compare("Parr1509_CP004010_SNPSummary.tsv", ["LOW", "MODIFIER", "MODERATE", "HIGH"], "binary", "kMeans")
    tSNEvsUMAP.plotCompare(df, "UMAP vs. tSNE mit k-Means-Clustering", "UMAP", "tSNE")
    tSNEvsUMAP.compareToMLST(df, 'cluster_allel_profiles.csv')
    print(sep)

    print("Mit SNP-Filtern nach MODERATE und HIGH, ohne MLST:")
    filter_df = tSNEvsUMAP.compare("Parr1509_CP004010_SNPSummary.tsv", ["MODERATE", "HIGH"], "binary", "kMeans")
    tSNEvsUMAP.plotCompare(filter_df, "UMAP vs. tSNE mit k-Means-Clustering (MODERATE, HIGH)", "UMAP", "tSNE")
    tSNEvsUMAP.compareToMLST(filter_df, 'cluster_allel_profiles.csv')
    print(sep)

    print("Ohne SNP-Filtern, mit MLST:")
    MLST_df = tSNEvsUMAP.compare("Parr1509_CP004010_SNPSummary.tsv", ["LOW", "MODIFIER", "MODERATE", "HIGH"], "binary", "kMeans",
                     ["TPANIC_RS00695", "TPANIC_RS02695", "TPANIC_RS03500"])
    tSNEvsUMAP.plotCompare(MLST_df, "MLST-UMAP vs. MLST-tSNE mit k-Means-Clustering", "MLST-UMAP", "MLST-tSNE")
    tSNEvsUMAP.compareToMLST(MLST_df, 'MLST_cluster_allel_profiles.csv')
    print(sep)

    print("Adjusted-Rand-Index der berechneten tSNEs und UMAPs:")
    tSNEvsUMAP.randIndex([df, filter_df, MLST_df])
    plt.show()