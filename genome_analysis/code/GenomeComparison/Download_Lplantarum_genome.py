import pandas as pd
import os


## configure
# wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
sp_id = "1590"    ## taxid of L.plantarum
levels = "Complete Genome"
to_download = "_genomic.fna.gz"
date = "20230208"

path2download = "../../data/"
path2download_dir = path2download + date + "_" + sp_id + "_fna"
if not os.path.exists(path2download_dir): os.mkdir(path2download_dir)

## set up data to download
data = pd.read_csv(path2download + date + "_assembly_summary.txt", sep="\t", skiprows=1)
data = data.loc[ (data["species_taxid"] == int(sp_id)) & (data["assembly_level"] == levels)  , :]
data.to_csv(path2download + date + "_info." + sp_id + "_genomic.csv")       ## save meta data

## download
for file_path in data["ftp_path"]:
    ac_name = file_path.split("/")[-1]
    file_name = "%s%s"  %(ac_name, to_download)
    os.system("wget " + file_path + "/" + file_name + " -P " + path2download_dir)
    os.system("gzip -d " + path2download_dir + "/" + file_name)



