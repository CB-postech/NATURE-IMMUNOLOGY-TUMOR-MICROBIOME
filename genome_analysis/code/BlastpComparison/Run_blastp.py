import pandas as pd
import os, sys

### setting path to blastdb
path2blastdb = "../../result/BlastpComparison/blastp/blastdb/"
if os.path.exists(path2blastdb): os.mkdir(path2blastdb)
path2blastdb = path2blastdb + "WCFS1"


### setting path to data & output
path2data = "../../data/BlastpComparison_1590-1622_faav2/"
db_strain = "Lp_WCFS1"
query_strain = ["ATCC8014", "ATCC14917", "ATCC35020", "IMB19", "Lp_DSM10667"]
path2output = path2blastdb.replace("blastdb/WCFS1", "")


### make blastdb
db_cmd = f"makeblastdb -in {path2data}{db_strain}_translated_cds.faa -parse_seqids -blastdb_version 5 -dbtype prot -out {path2blastdb}"
os.system(db_cmd)


### 
for i, query in enumerate(query_strain):
    path2input = path2data + query + "_translated_cds.faa"
    blast_cmd = f"blastp -db {path2blastdb} -query {path2input} -outfmt \"7 std qcovs qcovhsp\" -out {path2output}{query}_blastresult.txt"

    print(blast_cmd + "\n")
    os.system(blast_cmd)


#os.system("blastp -db ./VFDB_blast/VFDB_core_prot -query " + path2input + " -outfmt \"7 std qcovs qcovhsp\" -out " + path2output + "/RawBlastTabularResult.txt")

