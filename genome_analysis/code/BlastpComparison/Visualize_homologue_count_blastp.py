import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage
import sys
plt.rcParams["pdf.fonttype"] = 42

id_cut = 90
path2data = "../../result/BlastpComparison/Parse_hierarchical_blastp/"

data = pd.read_csv(path2data + "1590-1622_hierarchical_" + str(id_cut) + ".csv", index_col="strain")
data.index = [ strain.split("(")[0] for strain in data.index ]
cps_cluster = pd.read_csv(path2data + "wcfs1_cps_gene_locusID_faa_version.txt", sep="\t", index_col=0)
cps_cluster = cps_cluster.loc[ data.columns ,:]
data.columns = [ cps_cluster.loc[col, "old_locus_tag"] for col in cps_cluster.index ]


### dendrogram
plt.figure(figsize=(33,17))
plt.subplots_adjust(bottom=0.15)

Z = linkage(data, method="average")
D = dendrogram(Z, leaf_rotation=90., leaf_font_size=8., labels=data.index)
D["ivl"] = [ name.split("(")[0] for name in D["ivl"] ]

IMB19_idx = D["ivl"].index("IMB19")

plt.title("CDS similarity of CPS-related genes (Reference to WCFS1)")
plt.savefig(path2data + "1590-1622_hierarchical_" + str(id_cut) + ".png")
plt.clf()

### heatmap
cluster_color_map = {"1.0": "silver", "2.0":"orange", "3.0":"cyan", "4.0":"g"}
cluster_color = []
for cluster in cps_cluster["gene_cluster"]: cluster_color.append(cluster_color_map[str(cluster)])

ax = sns.clustermap(data, cmap="Reds", row_cluster=True, col_cluster=False, row_colors=None, col_colors=cluster_color,
                    method="average", figsize=(10,7), linewidths=.1, xticklabels=1, yticklabels=1)
plt.savefig(path2data + "1590-1622_hierarchicalHeatmap_" + str(id_cut) + ".png")
plt.savefig(path2data + "1590-1622_hierarchicalHeatmap_" + str(id_cut) + ".pdf", format="pdf")
plt.clf()



