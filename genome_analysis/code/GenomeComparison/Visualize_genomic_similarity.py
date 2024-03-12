import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram, cut_tree
import matplotlib.pyplot as plt


sp_id = "1590"
date = "20230208"
path2result = "../../result/" + date + "/"
path2data = "../../data/"
target_strains = ["VHProbi O04(GCF_024758665)", "ATCC 8014(GCF_002749655)", "TK-P2A(GCF_015377525)", "DOML(GCF_000604105)", "JDM1(GCF_000023085)", "BK-021(GCF_013487805)", "IMB19(IMB19)"]

genome_info = pd.read_csv(path2data + date + "_assembly_summary.txt", sep="\t", skiprows=1)
genome_info = genome_info.loc[ ( genome_info["species_taxid"] == int(sp_id) ) & ( genome_info["assembly_level"] == "Complete Genome" ),
                        ["# assembly_accession", "infraspecific_name"]]
genome_info["infraspecific_name"] = [ str(name).strip("strain=") for name in genome_info["infraspecific_name"] ]
genome_info["# assembly_accession"] = [ acc.split(".")[0] for acc in genome_info["# assembly_accession"] ]
genome_info.loc[ len(genome_info)+1 ] = ["IMB19", "IMB19"]
genome_info = genome_info.set_index("# assembly_accession").to_dict()["infraspecific_name"]
# genome_info = genome_info.set_index("infraspecific_name").to_dict()["# assembly_accession"]



data = pd.read_csv(path2result + "/dist/OrthoANI_lp.txt", sep="\t", index_col=0)
data.index = [ genome_info[acc] + "(" + acc + ")" if genome_info[acc] != "" else "unknown (" + acc + ")" for acc in data.index ]         ### need modification strain(ID)
data.columns = [ genome_info[acc] + "(" + acc + ")" if genome_info[acc] != "" else "unknown (" + acc + ")" for acc in data.columns ]     ### need modification strain(ID)

plt.figure(figsize=(30,15))
Z = linkage(data, method="average")
D = dendrogram(Z, leaf_rotation=90., leaf_font_size=9., labels=data.columns, distance_sort="descending")
IMB19_idx = D["ivl"].index("IMB19(IMB19)")

plt.subplots_adjust(bottom=0.20)
plt.title("genomic similarity of L. Plantarum \nIMB19: " + str(IMB19_idx))
plt.savefig(path2result + "/dist/LP_genomic_similarity.png")
plt.savefig(path2result + "/dist/LP_genomic_similarity.eps", format="eps")
plt.clf()


plt.figure(figsize=(7, 15))
data_subset = data.loc[target_strains, target_strains]
sub_Z = linkage(data_subset, method="average")
sub_D = dendrogram(sub_Z, leaf_rotation=90., leaf_font_size=9., labels=data_subset.columns, distance_sort="descending")

plt.subplots_adjust(bottom=0.15)
plt.title("Subset genomic similarity of L. Plantarum \nIMB19: " + str(IMB19_idx))
plt.yticks(np.arange(0.0001, 0.001, 0.0001))
plt.savefig(path2result + "/dist/LPsubset_genomic_similarity.png")
plt.savefig(path2result + "/dist/LPsubset_genomic_similarity.eps", format="eps")
plt.clf()
