import pandas as pd
import numpy as np


## configure
sp_id = "1590"
cdhit_cut = 50
cdhit_cuts = np.arange(50, 101, 5)
path2data = "../../data/"
path2wcfs1_data = path2data + "WCFS1/"
date = "20230208"
path2cdhit = "../../result/" + date + "/cdhit/"
path2result = "../../result/" + date + "/Parse_hierarchical/"

## find block: find a cluster of cdhit result containing cps gene locus ID
## contigs: dictionary of cdhit result
def find_block(gene, contigs):
    clutser_list = []
    for cluster in contigs.keys():
        if gene in contigs[cluster]: clutser_list.append( cluster )
    if len( clutser_list ) > 1: raise ValueError("one reference gene is located in different clusters")
    return clutser_list[0]

## load reference genes, refseq, locus_Tag
f = open(path2wcfs1_data + "WCFS1_CPSrelated_translated_cds.faa", "r")
data = f.read()
f.close()
blocks = data.split(">")
refs = {}
for block in blocks:
    refseq = block.split()[0]
    locus = [ tag.strip("[]").split("=")[-1] for tag in block.split() if tag.startswith("[locus_tag=") ]
    if len(locus) > 1: raise ValueError("refseqID and locus ID are not one to one match")
    refs[refseq] = locus

## load cluster ID and incorporate with refs seq
data = pd.read_csv(path2wcfs1_data + "wcfs1_cps_gene_locusID.txt", sep="\t")
refs = pd.DataFrame.from_dict(refs, orient="index", columns=["locus_tag"]).reset_index()
refs = pd.merge(left=refs, right=data, on="locus_tag", how="inner")
#refs.to_csv(path2result + "wcfs1_cps_gene_locusID_faa_version.txt", sep="\t", index=False)


## load strain info and faa ID info
strains = pd.read_csv(path2data + date + "_assembly_summary.txt", sep="\t", skiprows=1)
strains = strains.loc[ ( strains["species_taxid"] == int(sp_id) ) & ( strains["assembly_level"] == "Complete Genome" ),
                        ["# assembly_accession", "infraspecific_name", "ftp_path"]]
strains["infraspecific_name"] = [ str(name).strip("strain=") for name in strains["infraspecific_name"] ]
strains["strain_ID"] = [ name + "(" + acc + ")" for name, acc in zip(strains["infraspecific_name"], strains["# assembly_accession"]) ]
strains["faa_filename"] = [ filename.split("/")[-1] + "_translated_cds.faa" for filename in strains["ftp_path"] ]
strains = strains.drop(columns="ftp_path")
strains.loc[ len(strains) + 1 ] = ["IMB19", "IMB19", "IMB19(IMB19)", "IMB19.faa"]
strains = strains.loc[ strains["infraspecific_name"] != "WCFS1",:]

strain_faaID = {}
for strain, faa_file in zip(strains["strain_ID"], strains["faa_filename"]):
    f = open(path2data + "/" + date + "_" + sp_id + "_faa/" + faa_file, "r")
    faa = f.read()
    f.close()
    faa = faa.split(">")
    strain_faaID[strain] = [ ID.split(" ")[0] for ID in faa if ID != ""]


## table with (strain) * (reference cps gene)
for cut in cdhit_cuts:
    print(cut)
    ## load cdhit result
    f = open(path2cdhit + "db" + str(cut) + ".clstr")
    data = f.read()
    f.close()

    ## parsing cdhit result into dictionary -> (key)cluster: [ (value) homologous locus_ID ]
    lines = data.split("\n")
    contigs = {}
    for line in lines:
        if line.startswith(">Cluster "): 
            cluster = line.split(" ")[-1]
            contigs[ cluster ] = []
        else:
            contigs[ cluster ].append( line.split(", >")[-1].split("... ")[0] )

    ## count orthologs, iterating over the reference genes
    CNT_list, matched_list = [], []
    for gene in refs["index"]:
        cluster = find_block(gene, contigs)
        matched = pd.DataFrame([ ";".join(list(set(strain_faaID[strain]) & set(contigs[cluster]))) for strain in strain_faaID.keys() ], 
                    columns=["target_faa_ID"])
        matched["strain"] = strain_faaID.keys()
        matched["reference_faaID"] = [gene] * len(matched)
        matched_list.append(matched)

        count_hg = [ len(set(strain_faaID[strain]) & set(contigs[cluster])) for strain in strain_faaID.keys() ]
        count_hg_df = pd.DataFrame(count_hg, index=strain_faaID.keys(), columns=[gene])
        CNT_list.append( count_hg_df )
    matched_df = pd.concat(matched_list, axis=0)
    matched_df.to_csv(path2result + sp_id + "_matched_" + str(cut) + ".csv", index=False)
    CNT = pd.concat(CNT_list, axis=1)
    CNT.to_csv(path2result + sp_id + "_hierarchical_" + str(cut) + ".csv", index_label="strain")



