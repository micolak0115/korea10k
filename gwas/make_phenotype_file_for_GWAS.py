# %%
import math
import os

import pandas as pd

# %%
path_pheno = "/BiO/Research/GeroExpressome/Results/PhenoAge4/KGP_sbp_glucose_age_rdw_creat_bun_phenoage_calibrated.txt"
path_metadata = "/BiO/Access/kyungwhan1998/genome/gwas/Resources/Data/Total_Sex_Info_10K_Ver6_PlinkProcessing_Without_Age.txt"

# %%
df_pheno = pd.read_csv(path_pheno, sep=",")
df_pheno["phenoage_calibrated_advance"] = df_pheno["phenoage_calibrated"] - df_pheno["age"]
df_pheno["SampleID"] = df_pheno["Sample_ID_New"].apply(lambda x: "-".join(x.split("-")[:-1]))
df_pheno = df_pheno.drop_duplicates(subset=["SampleID"], keep="last")
df_pheno = df_pheno.drop(columns=["Sample_ID_New"])
df_pheno = df_pheno.set_index("SampleID")
list_col_pheno = list(df_pheno.columns)
df_pheno["#FID"] = df_pheno.index
df_pheno["IID"] = df_pheno.index
df_pheno = df_pheno[["#FID", "IID"] + list_col_pheno]
path_pca = "/BiO/Access/kyungwhan1998/genome/gwas/Resources/Data/9000Korean_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.excesshet_60.abhet_0.4.abhom_0.1.adsupport_0.9.kinship_3rd.prune.200kb_0.5.pca.eigenvec"
df_pca = pd.read_csv(path_pca, sep="\t")
list_samples_unrelated = df_pca["IID"]
list_samples_intsc = list(set(df_pheno.index).intersection(set(list_samples_unrelated)))
df_pheno_unrelated = df_pheno.loc[list_samples_intsc]
df_pheno_unrelated = df_pheno_unrelated.fillna("NA")

# %%
mean_val = df_pheno_unrelated["phenoage_calibrated_advance"].mean(axis=0)
stdev_val = df_pheno_unrelated["phenoage_calibrated_advance"].std(axis=0)
df_pheno_unrelated["phenoage_calibrated_advance_std"] = (df_pheno_unrelated["phenoage_calibrated_advance"] - mean_val) / stdev_val

# %%
median_val = df_pheno_unrelated["phenoage_calibrated_advance"].median(axis=0)
lower_quartile = np.percentile(df_pheno_unrelated["phenoage_calibrated_advance"], 25)
upper_quartile = np.percentile(df_pheno_unrelated["phenoage_calibrated_advance"], 75)
min_val = df_pheno_unrelated["phenoage_calibrated_advance"].min(axis=0)
max_val = df_pheno_unrelated["phenoage_calibrated_advance"].max(axis=0)
# %%
import matplotlib.pyplot as plt
import seaborn as sns

plt.figure(figsize=(5, 2))

sns.boxplot(
    data=df_pheno_unrelated,
    x="phenoage_calibrated_advance",
    width=0.15,
    color="lightgrey",
    linewidth=0.7,
    flierprops={"marker": "o", "markersize": 5, "alpha": 0.4}
)

plt.xlabel("PhenoAge Advance (calibrated)")
plt.ylabel("")
plt.tight_layout()
plt.show()

# %%
import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata


def inverse_normal_transform(x):
    x = pd.Series(x)
    mask = x.isna()
    ranks = rankdata(x[~mask], method='average')
    n = len(ranks)
    scaled_ranks = (ranks - 0.5) / n
    
    x_int = np.empty(len(x))
    x_int[mask] = np.nan
    x_int[~mask] = norm.ppf(scaled_ranks)
    
    return x_int
# %%

df_pheno_unrelated["phenoage_calibrated_advance_INT"] = inverse_normal_transform(df_pheno_unrelated["phenoage_calibrated_advance"])

# %%
df_pheno_unrelated.to_csv("/BiO/Access/kyungwhan1998/genome/gwas/Resources/Data/KGP_sbp_glucose_age_rdw_creat_bun_phenoage_calibrated_for_GWAS.txt", sep="\t", index=False)

# %%
# %%
