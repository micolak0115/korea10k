# %%
import os

import pandas as pd

# %%
path_clinical = "/BiO/Access/kyungwhan1998/genome/gwas/Resources/Data/KGP_sbp_glucose_age_rdw_creat_bun_phenoage_calibrated_for_GWAS.txt"
path_metadata = "/BiO/Access/kyungwhan1998/genome/gwas/Resources/Data/Total_Sex_Info_10K_Ver6_PlinkProcessing_Without_Age.txt"
path_pca = "/BiO/Access/kyungwhan1998/genome/gwas/Resources/Data/9000Korean_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.excesshet_60.abhet_0.4.abhom_0.1.adsupport_0.9.kinship_3rd.prune.200kb_0.5.pca.eigenvec"

# %%
def convert_sex_to_integer(sex):
    if sex == "M":
        sex = 1
    else:
        sex = 2
    
    return sex

# %%
df_clinical = pd.read_csv(path_clinical, sep="\t")[["#FID", "age"]]
df_clinical["age_sq"] = df_clinical["age"].apply(lambda x: x**2)
df_metadata = pd.read_csv(path_metadata, sep="\t", header=None, names=["#FID", "IID", "sex"])[["#FID", "sex"]]
df_metadata["sex"] = df_metadata["sex"].apply(convert_sex_to_integer)
df_pca = pd.read_csv(path_pca, sep="\t").iloc[:, :-10]
df_meta_pca = pd.merge(df_pca, df_metadata, how="inner", on=["#FID"])
df_meta_pca_clin = pd.merge(df_meta_pca, df_clinical, how="inner", on=["#FID"])
df_meta_pca_clin.to_csv("/BiO/Access/kyungwhan1998/genome/gwas/Resources/Data/KGP_age_sex_bmi_pc1_pc10_for_GWAS.txt", sep="\t", index=False)
df_meta_pca_clin

# %%
