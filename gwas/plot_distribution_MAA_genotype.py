# %%
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# %%
path_assoc = "/BiO/Access/kyungwhan1998/genome/gwas/Results/KGP_creat_age_glucose_sbp_phenoage_calibrated_GWAS.phenoage_calibrated_advance_INT.glm.linear"
gwas = pd.read_csv(path_assoc, sep="\t")
gwas_p_sorted = gwas[gwas["P"] < 1e-6].sort_values(by=["P"])
gwas_p_sorted["VAR"] = gwas_p_sorted["ID"] + "_" + gwas_p_sorted["REF"] + "_" + gwas_p_sorted["ALT"]

# %%
import subprocess

plink2 = "/BiO/Share/Tool/plink2"
dir_in = "/BiO/Access/kyungwhan1998/genome/pca/Resources/Data/ku10k"
dir_out = "/BiO/Access/kyungwhan1998/genome/gwas/Resources/Data"
bfile = os.path.join(dir_in, "9000Korean_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.excesshet_60.abhet_0.4.abhom_0.1.adsupport_0.9.kinship_3rd")
snplist = os.path.join(dir_out, "snplist.txt")
outfile = os.path.join(dir_out, "9000Korean_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.excesshet_60.abhet_0.4.abhom_0.1.adsupport_0.9.kinship_3rd.txt")

cmd = f"{plink2} \
  --bfile {bfile} \
  --extract {snplist} \
  --export A \
  --out {outfile}"
  
# subprocess.run(cmd, shell=True)

# %%
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# %%
path_genotype = "/BiO/Access/kyungwhan1998/genome/gwas/Resources/Data/9000Korean_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.excesshet_60.abhet_0.4.abhom_0.1.adsupport_0.9.kinship_3rd.txt.raw"
df_genotype = pd.read_csv(path_genotype, sep=r"\s+")

# %%
path_phenotype = "/BiO/Access/kyungwhan1998/genome/gwas/Resources/Data/KGP_creat_age_glucose_sbp_phenoage_calibrated_for_GWAS.txt"

df_phenotype = pd.read_csv(path_phenotype, sep="\t")[["IID", "phenoage_calibrated_advance_INT"]]

# %%
df_genotype_phenotype = df_genotype.merge(df_phenotype, how="inner", on="IID")

list_snps = list(filter(lambda x: str(x).startswith("chr"), list(df_genotype_phenotype.columns)))

import matplotlib.pyplot as plt
import numpy as np
# %%
import pandas as pd
import seaborn as sns
from scipy.stats import linregress

pheno_col = "phenoage_calibrated_advance_INT"

fig, axes = plt.subplots(1, 4, figsize=(10, 5), sharey=True)

for snp, ax in zip(list_snps, axes):
    d = df_genotype_phenotype[[pheno_col, snp]].replace([-9, -99], pd.NA).dropna()
    d[snp] = d[snp].astype(int)
    sns.violinplot(x=d[snp].astype("category"),
                   y=d[pheno_col],
                   data=d, 
                   inner=None, 
                   cut=0,
                   color="grey",
                   ax=ax)
    sns.boxplot(x=d[snp].astype("category"),
                y=d[pheno_col],
                data=d, 
                width=0.2,
                showcaps=True, 
                boxprops={
                        'zorder': 2, 
                        'color': "lightgrey",
                        },
                showfliers=False, 
                ax=ax)

    # ---------- Linear regression ----------
    slope, intercept, r_value, p_value, std_err = linregress(d[snp], d[pheno_col])

    # Points for regression line
    x_vals = np.array([0, 1, 2])
    y_vals = intercept + slope * x_vals
    ax.plot(x_vals, y_vals, linestyle="-", linewidth=5, color="firebrick")

    # Titles and labels
    ax.set_title(snp, fontsize=15)
    ax.set_xlabel("Genotype", fontsize=15)
    ax.set_ylabel("Mortality Age Acceleration\n(Inverse Normal Transformed)", fontsize=15)

plt.tight_layout()
plt.show()
plt.close()