import os
import sys
import glob
import palantir
import pandas as pd
import scanpy as sc
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager
import random
import pickle

%matplotlib inline

data_dir = "/home/dgcha/project/TIL/data/"
os.chdir(data_dir)
palantir_dir = os.path.expanduser(os.getcwd())

#import hvg_logcount and allgene_logcount
norm_df = pd.read_csv(os.path.join(data_dir, "logcd/til_logcd_hvg.csv"), sep=",", index_col=0)
norm_df_allgene = pd.read_csv(os.path.join(data_dir, "logcd/til_logcd.csv"), sep=",", index_col=0)

pcadim = 50
dm = 30

pca_projections, _ = palantir.utils.run_pca(norm_df, n_components=pcadim)
dm_res = palantir.utils.run_diffusion_maps(pca_projections, n_components=dm)
ms_data = palantir.utils.determine_multiscale_space(dm_res)

tsne = palantir.utils.run_tsne(ms_data, perplexity = 200)

fig, ax = palantir.plot.plot_tsne(tsne)

imp_df = palantir.utils.run_magic_imputation(norm_df_allgene, dm_res)

start_cell = np.argmax(imp_df['Ccr7'])
palantir.plot.highlight_cells_on_tsne(tsne, start_cell)
pr_res = palantir.core.run_palantir(ms_data, start_cell, num_waypoints=500)
palantir.plot.plot_palantir_results(pr_res, tsne)

#save
imp_df.to_csv('imp_df_til.csv')
ms_data.to_csv('ms_data_til.csv')
tsne.to_csv(data_dir + "/palantir_tsne_til.csv", index=True, index_label=True)
with open('pr_res_til.pickle', 'wb') as f:
    pickle.dump(pr_res, f, pickle.HIGHEST_PROTOCOL)

pr_res.branch_probs.to_csv("palantir_branch_probs_til.csv", index=True, index_label=True)
pr_res.pseudotime.to_csv("palantir_pseudotime_til.csv", index=True, index_label=True)