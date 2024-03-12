import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage
import sys

if len(sys.argv) > 1: cdhit_cut = sys.argv[1]
else:
    cdhit_cut = 50
    
target_strains = ["VHProbi O04(GCF_024758665)", "ATCC 8014(GCF_002749655)", "TK-P2A(GCF_015377525)", "DOML(GCF_000604105)", "JDM1(GCF_000023085)", "BK-021(GCF_013487805)", "IMB19(IMB19)"]
date = "20230208"
path2data = "../../result/" + date + "/Parse_hierarchical/"

data = pd.read_csv(path2data + "1590_hierarchical_" + str(cdhit_cut) + ".csv", index_col="strain")
data.index = [ strain[:-3] + ")" for strain in data.index ]
data.index = [ "unknown" + strain if strain.startswith("(") else strain for strain in data.index ]
data.index = [ strain.replace("(IMB)", "(IMB19)") if strain.startswith("IMB19") else strain for strain in data.index ]
cps_cluster = pd.read_csv(path2data + "wcfs1_cps_gene_locusID_faa_version.txt", sep="\t", index_col=0)
cps_cluster = cps_cluster.loc[ data.columns ,:]


### dendrogram
plt.figure(figsize=(30,15))
plt.subplots_adjust(bottom=0.20)

Z = linkage(data, method="average")
D = dendrogram(Z, leaf_rotation=90., leaf_font_size=9., labels=data.index, distance_sort="descending")
IMB19_idx = D["ivl"].index("IMB19(IMB19)")

plt.title("CDS similarity of CPS-related genes (Reference to WCFS1)\nIMB19: " + str(IMB19_idx))
plt.savefig(path2data + "1590_hierarchical_" + str(cdhit_cut) + ".png")
plt.savefig(path2data + "1590_hierarchical_" + str(cdhit_cut) + ".eps", format="eps")
plt.clf()


plt.figure(figsize=(10,12))
plt.subplots_adjust(bottom=0.20)

targets = D["ivl"][IMB19_idx-10:IMB19_idx+2]
data_subset = data.loc[ targets , :]
sub_Z = linkage(data_subset, method="average")
sub_D = dendrogram(sub_Z, leaf_rotation=90., leaf_font_size=9., labels=data_subset.index, distance_sort="descending")

plt.title("Subset: CDS similarity of CPS-related genes (Reference to WCFS1)\nIMB19: " + str(IMB19_idx))
plt.savefig(path2data + "1590_subsethierarchical_" + str(cdhit_cut) + ".png")
plt.savefig(path2data + "1590_subsethierarchical_" + str(cdhit_cut) + ".eps", format="eps")
plt.clf()





'''
### heatmap
cluster_color_map = {"1.0": "silver", "2.0":"orange", "3.0":"cyan", "4.0":"g"}
cluster_color = []
for cluster in cps_cluster["gene_cluster"]: cluster_color.append(cluster_color_map[str(cluster)])

ax = sns.clustermap(data, cmap="Reds", row_cluster=True, col_cluster=False, row_colors=None, col_colors=cluster_color,
                    method="average", figsize=(17,30), linewidths=.1, xticklabels=False, yticklabels=1)
plt.savefig(path2data + "1590_hierarchicalHeatmap_" + str(cdhit_cut) + ".png")
plt.clf()


data_focus = data.loc[ focus_strain ,:]
ax = sns.clustermap(data_focus, cmap="Reds", row_cluster=True, col_cluster=False, row_colors=None, col_colors=cluster_color,
                    method="average", figsize=(17,5), linewidths=.1, xticklabels=False, yticklabels=1)
plt.savefig(path2data + "1590_hierarchicalHeatmap_focus_strain_" + str(cdhit_cut) + ".png")
plt.clf()
'''

