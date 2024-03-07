import pandas as pd
import numpy as np
from scipy import stats
import os
import warnings
import matplotlib.pyplot as plt
from sklearn.model_selection import GridSearchCV, LeaveOneOut
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

# Ignore warnings
warnings.filterwarnings('ignore')

# Initialize constants
iternum = 100
test_size = 0.2
model = LogisticRegression()
param_grid = {'penalty': ['l2'], 'max_iter': [1e10], 'solver': ['lbfgs'],
              'C': np.arange(.1, 1., .1), 'class_weight': ['balanced']}

# Prepare data
scaler = StandardScaler()
edf = pd.read_csv("../data/EXP_IMVIGOR.txt", sep="\t")
edf = edf.set_index("NAME")
edf[edf.columns] = scaler.fit_transform(edf[edf.columns])
edf = edf.T

with open("../data/SigGenes.txt") as f:
    TGS_ = f.read().split()

tgs_common = set(edf.columns) & set(TGS_)
edf = edf[tgs_common]

TMB_col = 'FMOne mutation burden per MB'
RSP_col = 'Best Confirmed Overall Response'
pData_ = pd.read_csv('../data/Patient_IMVIGOR.txt')
pData_ = pData_.set_index("sample_id")
pData_ = pData_.assign(rsp=pd.NaT)
pData_.rsp.loc[pData_[RSP_col].isin(['PR', 'CR'])] = 1
pData_.rsp.loc[pData_[RSP_col].isin(['SD', 'PD'])] = 0
pData_ = pd.concat([pData_, edf], axis=1)
pData_ = pData_.dropna(subset=[TMB_col, 'rsp'] + list(tgs_common))

TGS = [tg for tg in TGS_ if tg in pData_.columns]
data = pData_.loc[pData_['Immune phenotype'] == 'inflamed']


# Prediction: TMB
X = data[[TMB_col]].to_numpy()
y = data.rsp.to_numpy().astype(np.int)
loo = LeaveOneOut()
gcv = GridSearchCV(model, param_grid=param_grid, cv=loo, scoring='accuracy', n_jobs=10)
gcv.fit(X, y)

prd = gcv.predict(X)
obs = y

tp, fp, tn, fn = 0, 0, 0, 0
for p, o in zip(prd, obs):
    if p == 1 and o == 1:
        tp += 1
    elif p == 1 and o == 0:
        fp += 1
    elif p == 0 and o == 1:
        fn += 1
    elif p == 0 and o == 0:
        tn += 1

_, p = stats.fisher_exact([[tp, fp], [fn, tn]])
OR = tp * tn / (fp * fn)
print("=== TMB only ===")
print("TP, TN, FP, FN: %d, %d, %d, %d"%(tp,tn,fp,fn))
print("OR: %f"%OR)
print("P: %f"%p)


# Prediction: TMB + SigGenes
n_tgs = 100 
tgs = TGS[:n_tgs]

X = data[[TMB_col] + list(set(data.columns) & set(tgs))].to_numpy()
y = data.rsp.to_numpy().astype(np.int)

loo = LeaveOneOut()
gcv = GridSearchCV(model, param_grid=param_grid, cv=loo, scoring='accuracy', n_jobs=10)
gcv.fit(X, y)

prd = gcv.predict(X)
obs = y

tp, fp, tn, fn = 0, 0, 0, 0
for p, o in zip(prd, obs):
    if p == 1 and o == 1:
        tp += 1
    elif p == 1 and o == 0:
        fp += 1
    elif p == 0 and o == 1:
        fn += 1
    elif p == 0 and o == 0:
        tn += 1

_, p = stats.fisher_exact([[tp, fp], [fn, tn]])
OR = tp * tn / (fp * fn)
print("=== TMB + SigGenes ===")
print("TP, TN, FP, FN: %d, %d, %d, %d"%(tp,tn,fp,fn))
print("OR: %f"%OR)
print("P: %f"%p)
