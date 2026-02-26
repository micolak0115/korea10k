# %%
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2

# %%
path_assoc = "/BiO/Access/kyungwhan1998/genome/gwas/Results/KGP_creat_age_glucose_sbp_phenoage_calibrated_GWAS.phenoage_calibrated_advance_INT.glm.linear"
path_vep = "/BiO/Research/Korea10KGenome/Results/VariantAnnotation.by.VEP115/Korea10K_VariantAnnotation.FilterIndividualABHet2STD.VariantABHet0.4ABHom0.9ADSupp0.9.Biallelic.Geno_0.01.Mind_0.1.HWE_1e6.Het_3STD.Kinship_3rd.nonKorean.20251114.maf"
# path_vep_gwas = "/BiO/Access/kyungwhan1998/genome/gwas/Results.phenoage_calibrated_advance.glm.linear.annot"
gwas = pd.read_csv(path_assoc, sep="\t")
# %%
gwas_p_sorted = gwas[gwas["P"] < 1e-6].sort_values(by=["P"])

# %%
gwas_p_sorted["VAR"] = gwas_p_sorted["ID"] + "_" + gwas_p_sorted["REF"] + "_" + gwas_p_sorted["ALT"]

# %%
def calc_lambda_median_from_df(df):
    """
    Compute genomic inflation factor lambda (median) from PLINK output DataFrame.
    We expect columns: BETA and SE or Z (ZSTAT). We'll use (BETA/SE)**2 if possible,
    otherwise try ZSTAT**2, otherwise compute from p-value (chi2 = qchisq(1-p,1)).
    """
    if 'BETA' in df.columns and 'SE' in df.columns:
        z = df['BETA'] / df['SE']
        chi2_vals = z**2
    elif 'ZSTAT' in df.columns:
        chi2_vals = df['ZSTAT']**2
    elif 'T_STAT' in df.columns:
        chi2_vals = df['T_STAT']**2
    elif 'P' in df.columns:
        # convert p to chi2: chi2 = quantile(1-p, df=1)
        p = df['P'].astype(float)
        # avoid p==0 (set small)
        p = p.clip(lower=1e-300)
        chi2_vals = chi2.isf(p, df=1)
    elif 'p' in df.columns:
        p = df['p'].astype(float)
        p = p.clip(lower=1e-300)
        chi2_vals = chi2.isf(p, df=1)
    else:
        raise ValueError("No suitable fields (BETA/SE or ZSTAT or P/p) found to compute lambda.")
    # Remove NaN/inf
    chi2_vals = chi2_vals.replace([np.inf, -np.inf], np.nan).dropna()
    if len(chi2_vals) == 0:
        return np.nan
    lambda_median = np.median(chi2_vals) / 0.456  # median chi2(1df) == 0.456
    return lambda_median

# %%
def qq_plot(df, lambda_median, ax=None):
    import numpy as np
    pvals = df["P"]
    if ax is None:
        fig, ax = plt.subplots(figsize=(3, 3))
    pvals = np.array(pvals.dropna().astype(float))
    pvals = pvals[(pvals>0) & (pvals<=1)]
    n = len(pvals)
    expected = -np.log10((np.arange(1, n+1) - 0.5) / n)
    observed = -np.log10(np.sort(pvals))
    ax.scatter(expected, observed, s=5, color="k")
    maxv = max(expected.max(), observed.max()) + 0.5
    ax.plot([0, maxv], [0, maxv], linestyle='--', color="firebrick")
    ax.text(0, 16, s=round(float(lambda_median), 3))
    ax.set_xlabel("Expected -log10(p)")
    ax.set_ylabel("Observed -log10(p)")
    return ax

# %%
lambda_median = calc_lambda_median_from_df(gwas)

# %%
import statsmodels.api as sm
from scipy import stats

fig, ax = plt.subplots(1, 1, figsize=(5, 5))
results = stats.probplot(gwas["P"], dist=stats.uniform, plot=None)
ax.scatter(-np.log10(results[0][0]), -np.log10(results[0][1]), s=5, color="k")
ax.plot([0, 8], [0, 8], color="firebrick")
ax.text(0, 8, s=f"$\\lambda_{{GC}}$ = {round(float(lambda_median), 3)}", fontsize=12)
ax.set_xlabel("Expected -log10(P)", fontsize=15)
ax.set_ylabel("Observed -log10(P)", fontsize=15)
plt.show()
plt.close()

