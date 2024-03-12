import os
import numpy as np
import pandas as pd

def word_size(t):
    if t >= 0.7:
        return 5
    elif t >= 0.6:
        return 4
    elif t >= 0.5:
        return 3
    elif t >= 0.4:
        return 2
    
    return False


## configure
path2cdhit = "/path/to/cdhit"
path2input = "../../data/20230208_1590_all-faa/all.faa"
path2output = "../../result/20230208/cdhit/db"
threshold_range = np.arange(100, 39, -1)


for threshold in threshold_range:
    c = "%.2f"%(float(threshold)/100)
    n = word_size(float(c))

    ## setting command
    cmd = path2cdhit
    cmd += " -i " + path2input
    cmd += " -o " + path2output + str(threshold)
    cmd += " -c " + c
    cmd += " -n " + str(n)
    cmd += " -aS 0.9" 
    cmd += " -M 32000"
    cmd += " -d 0"
    cmd += " -T 16"

    ## run cd-hit
    #print(cmd)
    os.system(cmd)






