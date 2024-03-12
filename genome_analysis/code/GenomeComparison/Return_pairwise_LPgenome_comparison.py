import pandas as pd
import numpy as np 
from itertools import combinations
import glob
import os, time, sys


## configure
path2result = "../../result/20230208/ani_Lplantarum/"
path2data = "../../data/20230208_1590_fna/"
path2OAT = "/path/to/OrthoANI/"
path2blast = "/path/to/blast"
bin_number = 15
print(sys.argv)
if len(sys.argv) == 1: current_bin = ""
if len(sys.argv) > 1 and sys.argv[1] == "-cb": current_bin = sys.argv[2]


## pair generation and divide bin
targets = glob.glob(path2data + "*.fna")
pairs = list(combinations(targets, 2))
pairs_df = pd.DataFrame([pairs], columns=pd.cut(np.arange(0, len(pairs)), bin_number, labels=False), index=["pairs"]).T
if not current_bin == "": pairs_df = pairs_df.loc[ int(current_bin) ,:]
print(pairs_df)


## run pairwise genome similarity calculation
idx = 1
for filei, filej in pairs_df["pairs"].tolist():
    tagi = filei.split("/")[-1].split(".")[0]    
    tagj = filej.split("/")[-1].split(".")[0]
    
    tag = tagi + "-" + tagj
    print(tag, time.ctime())

    cmd = "java -jar " + path2OAT + " -num_threads 4 -blastplus_dir " + path2blast + " -fasta1 " + filei + " -fasta2 " + filej + " > " + path2result + tag + ".txt"
    #print(cmd)
    os.system(cmd)

    ## check process
    print( str( (idx / len(pairs))*100 ) + "% (" + str(idx) + "/" + str(len(pairs_df)) + ")" )
    idx += 1

print("all pairs: " + str(len(pairs)))

