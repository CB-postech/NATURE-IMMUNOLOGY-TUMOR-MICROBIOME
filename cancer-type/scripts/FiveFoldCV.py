import pandas as pd
import numpy as np
from scipy import stats
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve,auc 
import warnings

warnings.filterwarnings('ignore')
warnings.filterwarnings(action='ignore', category=DeprecationWarning)
warnings.filterwarnings(action='ignore', category=FutureWarning)

# Initialize
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
pfs_tmb = []

for idx in range(iternum):
    X = data[[TMB_col]].to_numpy()
    y = data.rsp.to_numpy().astype(np.int)

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=idx, stratify=y)

    gcv = GridSearchCV(model, param_grid=param_grid, cv=5, scoring='roc_auc', n_jobs=10).fit(X_train, y_train)
    prd = gcv.best_estimator_.predict(X_test)
    fpr, tpr, threshold = roc_curve(y_test, prd, pos_label=1)
    AUC = auc(fpr, tpr)
    pfs_tmb.append(AUC)

# Predict: +mac
pfs = []

tgs = TGS[:100]
for idx in range(iternum):
    X = data[[TMB_col] + list(set(data.columns) & set(tgs))].to_numpy()
    y = data.rsp.to_numpy().astype(np.int)

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=idx, stratify=y)

    gcv = GridSearchCV(model, param_grid=param_grid, cv=5, scoring='roc_auc', n_jobs=10).fit(X_train, y_train)
    prd = gcv.best_estimator_.predict(X_test)
    fpr, tpr, threshold = roc_curve(y_test, prd, pos_label=1)
    AUC = auc(fpr, tpr)
    pfs.append(AUC)


print("AUC, TMB only: %f"%(stats.tmean(pfs_tmb)))
print("AUC, TMB + SigGenes: %f"%(stats.tmean(pfs)))
print("P: %f"%stats.ttest_ind(pfs_tmb, pfs)[1])

fw = open("../results/AUCs.txt", 'w')
for a, b in zip(pfs_tmb, pfs):
    fw.write("%s\t%s\n" % (a, b))
fw.close()

