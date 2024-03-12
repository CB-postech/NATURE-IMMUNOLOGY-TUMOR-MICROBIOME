import pandas as pd
import numpy as np
import os, sys

sp_id = "1590-1622"
cut = 90
query_strain = ["ATCC8014", "ATCC14917", "ATCC35020", "IMB19", "Lp_DSM10667"]
path2output = "/home/jhlee/ProjectArchives/genome_comparison_Lplantarum_strain/result/BlastpComparison/blastp/"
path2result = path2output.replace("blastp", "Parse_hierarchical_blastp")

wcfs1_cps = pd.read_csv("../../result/BlastpComparison/Parse_hierarchical/wcfs1_cps_gene_locusID_faa_version.txt",
                        sep="\t")


h_result_list = []
for i, query in enumerate(query_strain):

    ## filter homologous proteins using specific cut-off
    result = pd.read_csv(path2output + f"{query}_blastresult.txt", sep="\t", comment="#", header=None)
    col = pd.read_csv(path2output + f"{query}_blastresult.txt", sep="\t", nrows=4, header=None).iloc[3].tolist()[0]
    col = col.split(": ")[-1].split(", ")
    result.columns = col
    result = result.loc[ result["% query coverage per hsp"] > 60 ,:]
    result = result.loc[ result["% identity"] > cut ]            ### parameter used by virulencefinder
    result = result.sort_values(by=["% identity", "bit score"], ascending=False)

    ## count homologous proteins
    print(query)
    homologues = []
    for pid_ in wcfs1_cps["index"]:
        pid = pid_.split("|")[-1]
        result_ = result.loc[result["subject acc.ver"] == pid, :]
        h_prot = pd.DataFrame([str(len(result_))], index=[pid_], columns=[query])
        print(pid + ": " + str(len(result_)))

        #if len(result_) == 0:
        #    matched_prot = pd.DataFrame([pid_, ""], index=["targetID", "subjectID"]).T
        #else:
        #    matched_prot = pd.DataFrame([pid_, ":".join(result_["query acc.ver"].tolist())], 
        #                                index=["targetID", "subjectID"]).T

        homologues.append(h_prot)
    homology_df = pd.concat(homologues, axis=0)

    h_result_list.append( homology_df )
h_result = pd.concat(h_result_list, axis=1)
h_result = h_result.T
h_result.to_csv(path2result + sp_id + "_hierarchical_" + str(cut) + ".csv", index_label="strain")