import pandas as pd
import numpy as np
import glob

## configure
path2result = "../../result/20230208/dist/"
path2ANI = "../../result/20230208/ani_Lplantarum/"
target_list = glob.glob(path2ANI + "*.txt")

sim_matrix_list = []
for target in target_list:
    f = open(target)
    data = f.read()
    f.close()

    data = data.split()
    try:
        data = data[ data.index("OrthoANI"): ]
        p = float(data[2])
    except:
        print(target)
        continue

    a,b = target.split("/")[-1][:-4].split("-")
        
    tmp_df = pd.DataFrame([[a, b, p/100], [b, a, p/100]], columns=["acc1", "acc2", "similarity"])
    sim_matrix_list.append(tmp_df)
sim_matrix = pd.concat(sim_matrix_list, axis=0).sort_values(by=["acc1", "acc2"])
mat = sim_matrix.pivot(index="acc1", columns="acc2", values="similarity")
mat = mat.fillna(1)
mat.to_csv(path2result + "OrthoANI_lp.txt", sep="\t")


