# %%
import itertools

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from statsmodels.stats.outliers_influence import variance_inflation_factor

# %%
path_pheno = "/BiO/Access/kyungwhan1998/genome/gwas/Resources/Data/KGP_creat_age_glucose_sbp_phenoage_calibrated_for_GWAS.txt"
path_covar = "/BiO/Access/kyungwhan1998/genome/gwas/Resources/Data/KGP_age_sex_bmi_pc1_pc10_for_GWAS.txt"

df_pheno = pd.read_csv(path_pheno, sep="\t").drop(columns=["IID"])
df_covar = pd.read_csv(path_covar, sep="\t").drop(columns=["IID", "age"])
df_pheno_covar = pd.merge(df_pheno, df_covar, how='inner', on="#FID")

# %%
import matplotlib.pyplot as plt
import seaborn as sns

plt.figure(figsize=(2, 4))
sns.boxplot(data=df_pheno_covar, x="sex", y="phenoage_calibrated_advance_INT")
plt.ylim(top=20)
plt.axhline(y=0, color="red", linewidth=2)
plt.xticks([0, 1], ["Male", "Female"])
plt.show()
plt.close()
# %%
# %%
# df_pheno_covar = pd.read_csv("/BiO/Access/kyungwhan1998/genome/gwas/Resources/Data/KGP_creat_age_glucose_sbp_phenoage_calibrated_for_GWAS.txt", sep="\t")
# cov_cols = ["sex", "age"] + [f"PC{i}" for i in range(1, 11)]
cols = list(df_pheno_covar.columns)[2:]
X = df_pheno_covar[cols].dropna().astype(float)

list_comb_cov = list(itertools.combinations(cols, r=2))

for c1, c2 in list_comb_cov:
    if c1 == "phenoage_calibrated_advance_INT" or c2 == "phenoage_calibrated_advance_INT":
        r, p = pearsonr(X[c1], X[c2])
        
        if p < 0.05:
            print(f"{c1} vs {c2}: r={r:.4f}, p={p:.4e}")

# %%
def compute_vif(df):
    return pd.Series(
        [variance_inflation_factor(df.values, i) for i in range(df.shape[1])],
        index=df.columns
    )

vif = compute_vif(X)
print("Initial VIF:\n", vif)

# Remove covariates with VIF > 50 (UKB threshold)
while vif.max() > 50:
    drop_col = vif.idxmax()
    print(f"Dropping collinear covariate: {drop_col}")
    X = X.drop(columns=[drop_col])
    vif = compute_vif(X)

print("\nFinal covariates to use:")
print(list(X.columns))