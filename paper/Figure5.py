# %%
import math
import os
import pickle

import lifelines
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import mortality
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
from adjustText import adjust_text
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy import stats
from scipy.stats import chi2
from statannotations.Annotator import Annotator


# %%
def load_pickle(pickle_file):
    with open(file=pickle_file, mode='rb') as f:
        obj = pickle.load(f)
    
    return obj

def plot_partial_hazard_into(ax, df, event_col='status', col_ph="partial_hazard"):
    g = sns.boxplot(
        data=df, 
        x=event_col, 
        y=col_ph, 
        color="grey",
        width=0.5,
        ax=ax
    )

    # --- Add annotation using statannotations ---
    pairs = [(0, 1)]
    annotator = Annotator(ax, pairs, data=df, x=event_col, y=col_ph)
    annotator.configure(test='Mann-Whitney', text_format='star', loc='inside', verbose=0)
    annotator.apply_and_annotate()

    # --- Labeling ---
    g.set_ylabel("$log_{10}$(Partial Hazard)", fontsize=15)
    g.set_xlabel("Death Status", fontsize=15)
    g.set_xticks([0, 1], ["Alive", "Deceased"])
    
def get_alive_dead_counts(df_partial, original_processed=None, status_col="status"):
    if status_col in df_partial.columns:
        counts = df_partial[status_col].value_counts(dropna=True)
        n_dead = int(counts.get(1, 0))
        n_alive = int(counts.get(0, 0))
        return n_alive, n_dead, n_alive + n_dead

    if original_processed is not None and status_col in original_processed.columns:
        intersect = df_partial.index.intersection(original_processed.index)
        if len(intersect) == 0:
            return len(df_partial), 0, len(df_partial)
        counts = original_processed.loc[intersect, status_col].value_counts(dropna=True)
        n_dead = int(counts.get(1, 0))
        n_alive = int(counts.get(0, 0))
        return n_alive, n_dead, n_alive + n_dead

    return len(df_partial), 0, len(df_partial)

def normalize_variant(var):
    ref = var.split("_")[-2]
    alt = var.split("_")[-1]
    other = "_".join(var.split("_")[:-2])
    
    if len(ref) > len(alt):
        print("deletion")
        ref_new = ref[1:]
        alt_new = "-"
        var = other + "_" + ref_new + "_" + alt_new
        print(var)
    
    if len(ref) < len(alt):
        print("insertion")
        ref_new = "-"
        alt_new = alt[1:]
        var = other + "_" + ref_new + "_" + alt_new
        print(var)
        
    return var

def manhattan_plot(
    df,
    chr_col='#CHROM',
    pos_col='POS',
    p_col='P',
    ax=None,
    sigline=5e-8,
    vep_df=None
    ):
    
    list_annot_text = []
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(12,4))

    df = df.copy()
    df = df[df[chr_col].notnull()]

    df[chr_col] = df[chr_col].astype(str).str.replace("chr", "", regex=False)

    chr_map = {'X': 23, 'Y': 24}
    df['CHRint'] = df[chr_col].replace(chr_map).astype(int)
    df['CHRstr'] = df['CHRint'].astype(str)

    lead_idx = df.groupby('CHRint')[p_col].idxmin()
    df['is_lead'] = False
    df.loc[lead_idx, 'is_lead'] = True

    positions = []
    ticks = []
    tick_labels = []
    cumul = 0

    for chr_id, group in df.groupby('CHRint'):
        group = group.sort_values(pos_col)
        group['cumuBP'] = group[pos_col] + cumul
        positions.append(group)

        mid = (group['cumuBP'].min() + group['cumuBP'].max()) / 2
        ticks.append(mid)
        tick_labels.append(str(chr_id))

        cumul = group['cumuBP'].max() + 1e6

    plot_df = pd.concat(positions)

    unique_chrs = sorted(plot_df['CHRstr'].unique(), key=int)
    colors = ["0.3" if i%2==0 else "0.6" for i in range(len(unique_chrs))]
    color_map = {c: colors[i] for i, c in enumerate(unique_chrs)}
    plot_df['color'] = plot_df['CHRstr'].map(color_map)

    plot_df['-log10p'] = -np.log10(plot_df[p_col].astype(float).clip(lower=1e-300))
    ax.scatter(plot_df['cumuBP'], plot_df['-log10p'], c=plot_df['color'], s=5)

    # -------------------------
    # Only annotate & plot lead SNPs
    # -------------------------
    lead_df = plot_df[plot_df['is_lead']]
    lead_df_suggestive = lead_df[lead_df[p_col] < sigline]

    ax.scatter(
        lead_df_suggestive['cumuBP'],
        lead_df_suggestive['-log10p'],
        c='firebrick',
        s=20
    )

    # Axes
    ax.set_xlabel("Chromosome")
    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels, rotation=90)
    ax.set_ylabel("-log10(P)")

    if sigline:
        ax.axhline(-np.log10(sigline), color='firebrick', linestyle='--', linewidth=1)

    # -------------------------------------------------------------
    # Annotation ONLY for suggestive lead SNPs
    # -------------------------------------------------------------
    if vep_df is not None:

        v = vep_df.copy()

        # Normalize
        v['Chromosome'] = (
            v['Chromosome'].astype(str)
            .str.replace("chr", "", regex=False)
            .replace({'X': 23, 'Y': 24})
        )

        v = v.dropna(subset=['Chromosome', 'vcf_pos'])
        v['Chromosome'] = v['Chromosome'].astype(int)
        v['vcf_pos'] = v['vcf_pos'].astype(int)

        # Merge ONLY suggestive lead SNPs
        merged = pd.merge(
            lead_df_suggestive,
            v,
            left_on=['CHRint', pos_col],
            right_on=['Chromosome', 'vcf_pos'],
            how='left'
        )

        ann=merged[merged['Annotation'].notna()]

        # Debug unmatched
        unmatched = merged[merged['Annotation'].isna()]
        if len(unmatched) > 0:
            print("Unmatched suggestive lead SNPs:")
            print(unmatched[[chr_col, pos_col, p_col, 'CHRint']].head())

        for _, row in ann.iterrows():
            annot_text = ax.text(
                row['cumuBP'],
                row['-log10p']+0.3,
                row['Annotation'],
                fontsize=plt.rcParams["font.size"],
                ha='center',
                color='firebrick',
                style='italic',
                weight='bold'
            )
            list_annot_text.append(annot_text)

    return ax, list_annot_text

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
    lambda_median=np.median(chi2_vals) / 0.456  # median chi2(1df) == 0.456
    return lambda_median

def qq_plot(df, lambda_median, ax=None):
    import numpy as np
    pvals = df["P"]
    if ax is None:
        fig, ax = plt.subplots(figsize=(3, 3))
    pvals = np.array(pvals.dropna().astype(float))
    pvals = pvals[(pvals>0) & (pvals<=1)]
    n=len(pvals)
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
plt.rcParams["font.size"] = 16
fig = plt.figure(figsize=(20, 15))

gs = gridspec.GridSpec(
    2, 4,
    width_ratios=[1, 1, 1, 1],
    height_ratios=[1, 0.5],
    wspace=0.5,
    hspace=0.35
)

ax_0 = fig.add_subplot(gs[0, 0])
ax_1 = fig.add_subplot(gs[0, 1])
ax_2 = fig.add_subplot(gs[0, 2])
ax_3 = fig.add_subplot(gs[0, 3])
ax_4 = fig.add_subplot(gs[1, 0:3])
ax_5 = fig.add_subplot(gs[1, 3])

labels = ["A", "B", "C", "D", "E", "F"]
axes = {
        ax_0:[-0.3, 1.1], 
        ax_1:[-0.3, 1.1], 
        ax_2:[-0.3, 1.1], 
        ax_3:[-0.3, 1.1],
        ax_4:[-0.075, 1.15],
        ax_5:[-0.3, 1.15]
        }

for idx, (ax, xy) in enumerate(axes.items()):
    ax.text(
        xy[0], xy[1], labels[idx],
        transform=ax.transAxes,
        fontsize=plt.rcParams["font.size"] + 10,
        fontweight="bold",
        va="top",
        ha="right"
    )

# --------------------------------------------------------------------
# (A) Forest plot
# --------------------------------------------------------------------
workdir = "/BiO/Research/GeroExpressome/Results/Cox_PH_Model"
path_dict_params = os.path.join(workdir, "dict_params_creat_waist_bmi_trig_age_glucose_sbp_hdl_lymph_dbp_totchol_monopa.pk")
path_cox_model = os.path.join(workdir, "cox_model_creat_waist_bmi_trig_age_glucose_sbp_hdl_lymph_dbp_totchol_monopa_min_plus_std_deviance.pk")

dict_params = load_pickle(path_dict_params)
cox_model = load_pickle(path_cox_model)
dict_cox_hr = dict(cox_model.hazard_ratios_)
dict_cox_hr_sig = {k: v for k, v in dict_cox_hr.items() if round(abs(1 - v), 4) > 0}
dict_cox_lower_ci = dict(cox_model.confidence_intervals_["95% lower-bound"])
dict_cox_exp_lower_ci = {k:np.exp(v) for k,v in dict_cox_lower_ci.items()}
dict_cox_exp_lower_ci_sig = {k:v for k,v in dict_cox_exp_lower_ci.items() if k in list(dict_cox_hr_sig)}
dict_cox_higher_ci = dict(cox_model.confidence_intervals_["95% upper-bound"])
dict_cox_exp_higher_ci = {k:np.exp(v) for k,v in dict_cox_higher_ci.items()}
dict_cox_exp_higher_ci_sig = {k:v for k,v in dict_cox_exp_higher_ci.items() if k in list(dict_cox_hr_sig)}
variables = list(dict_cox_hr.keys())
data = sorted(
    [(v, dict_cox_hr[v], dict_cox_exp_lower_ci[v], dict_cox_exp_higher_ci[v]) for v in variables],
    key=lambda x: x[1]
)

labels = [d[0] for d in data]
hr_vals = [d[1] for d in data]
ci_low = [d[2] for d in data]
ci_high = [d[3] for d in data]

y = np.arange(len(labels))
ax_0.set_yticks(y)
ax_0.set_yticklabels(labels, fontsize=plt.rcParams["font.size"]+6)

for idx, (label, hr, low, high) in enumerate(data):
    color = "firebrick" if round(float(hr), 4) > 1.0000 else "black"
    ax_0.hlines(idx, low, high, color=color, linewidth=3)
    ax_0.plot(hr, idx, "d", markersize=10, color=color, zorder=3)
    if round(float(hr), 4) > 1.0000:
        ax_0.text(high + 0.02, idx, round(float(hr), 3),
                  color="firebrick", va="center", fontsize=plt.rcParams["font.size"]+2, weight="bold")

for tick_label, hr in zip(ax_0.get_yticklabels(), hr_vals):
    if round(float(hr), 4) > 1.000:
        tick_label.set_color("firebrick")
        tick_label.set_fontweight("bold")

c_index = c_index = round(float(cox_model.concordance_index_), 3)
ax_0.text(1.4, 0, s=f"C-index\n= {c_index}", fontsize=plt.rcParams["font.size"]+4, weight="bold")
ax_0.axvline(1, linestyle="--", color="grey", linewidth=1)
ax_0.set_xlabel("Hazard Ratio (HR)", fontsize=plt.rcParams["font.size"]+6)
ax_0.set_xticks(np.arange(0.8, 2.0, 0.2))
sns.despine(ax=ax_0)
ax_0.grid(axis="x", linestyle="-", alpha=0.5)

# --------------------------------------------------------------------
# (B) NHANES III
# --------------------------------------------------------------------
path_nhanes_iii = "/BiO/Research/GeroExpressome/Resources/Data/NHANES/NHANES3.txt"
path_nhanes_iv = "/BiO/Research/GeroExpressome/Resources/Data/NHANES/NHANES4.txt"
path_kgp = f"/BiO/Research/GeroExpressome/Results/PreprocessedData/KGP_biomarkers_only_common_NHANES3.txt"

outfilename = "log_transformed_partial_hazard_{0}_{1}.tsv"

nhanes_iii_data = pd.read_csv(path_nhanes_iii)
nhanes_iii_processed_data = mortality.parse_nhanes_data_mortality_stats(nhanes_iii_data)
nhanes_iii_processed_data = nhanes_iii_processed_data.rename(columns={"sampleID": "Sample_ID_New"})
nhanes_iii_processed_data = nhanes_iii_processed_data.set_index("Sample_ID_New")

nhanes_iv_data = pd.read_csv(path_nhanes_iv)
nhanes_iv_processed_data = mortality.parse_nhanes_data_mortality_stats(nhanes_iv_data)
nhanes_iv_processed_data = nhanes_iv_processed_data.rename(columns={"sampleID": "Sample_ID_New"})
nhanes_iv_processed_data = nhanes_iv_processed_data.set_index("Sample_ID_New")

kgp_data = pd.read_csv(path_kgp, sep="\t")
kgp_data = kgp_data.set_index("Sample_ID_New")

dict_biomarkers_hazard = mortality.get_biomarkers_hazard(cox_model)
biomarkers = list(dict_biomarkers_hazard.keys())

df_partial_hazard_nhanes_iii = mortality.calculate_partial_hazard(
    cox_model, nhanes_iii_processed_data.reset_index(), biomarkers, col_status="status"
).set_index("Sample_ID_New")

df_partial_hazard_nhanes_iv = mortality.calculate_partial_hazard(
    cox_model, nhanes_iv_processed_data.reset_index(), biomarkers, col_status="status"
).set_index("Sample_ID_New")

df_partial_hazard_kgp = mortality.calculate_partial_hazard(
    cox_model, kgp_data.reset_index(), biomarkers
).set_index("Sample_ID_New")

n_alive_iii, n_dead_iii, n_total_iii = get_alive_dead_counts(
    df_partial_hazard_nhanes_iii, nhanes_iii_processed_data
)
n_alive_iv, n_dead_iv, n_total_iv = get_alive_dead_counts(
    df_partial_hazard_nhanes_iv, nhanes_iv_processed_data
)


plot_partial_hazard_into(ax_1, df_partial_hazard_nhanes_iii, col_ph="log10(partial_hazard)")
ax_1.set_title(f"NHANES III\n(N={len(df_partial_hazard_nhanes_iii):,})",
               fontsize=plt.rcParams["font.size"]+6)
ax_1.set_ylabel(ax_1.get_ylabel(), fontsize=plt.rcParams["font.size"]+6)
ax_1.set_xlabel(ax_1.get_xlabel(), fontsize=plt.rcParams["font.size"]+6)
ax_1.set_xticks([0, 1])
ax_1.set_xticklabels([f"Alive\n(N={n_alive_iii:,})", f"Deceased\n(N={n_dead_iii:,})"])
ax_1.grid(axis="y", linestyle="-", alpha=0.5)
sns.despine(ax=ax_1)

# --------------------------------------------------------------------
# (C) NHANES IV
# --------------------------------------------------------------------

plot_partial_hazard_into(ax_2, df_partial_hazard_nhanes_iv, col_ph="log10(partial_hazard)")
ax_2.set_title(f"NHANES IV\n(N={len(df_partial_hazard_nhanes_iv):,})",
               fontsize=plt.rcParams["font.size"]+6)
ax_2.set_ylabel(ax_2.get_ylabel(), fontsize=plt.rcParams["font.size"]+6)
ax_2.set_xlabel(ax_2.get_xlabel(), fontsize=plt.rcParams["font.size"]+6)
ax_2.set_xticks([0, 1])
ax_2.set_xticklabels([f"Alive\n(N={n_alive_iv:,})", f"Deceased\n(N={n_dead_iv:,})"])
ax_2.grid(axis="y", linestyle="-", alpha=0.5)
sns.despine(ax=ax_2)

# --------------------------------------------------------------------
# (D) KGP PAA Boxplot
# --------------------------------------------------------------------
path_pheno = "/BiO/Research/GeroExpressome/Results/PhenoAge4/KGP_creat_age_glucose_sbp_phenoage_calibrated.txt"
path_metadata = "/BiO/Access/kyungwhan1998/genome/gwas/Resources/Data/Total_Sex_Info_10K_Ver6_PlinkProcessing_Without_Age.txt"

df_meta = pd.read_csv(path_metadata, sep="\t", header=None, names=["#FID", "IID", "Sex"])[["IID", "Sex"]]

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

mean_val = df_pheno_unrelated["phenoage_calibrated_advance"].mean(axis=0)
stdev_val = df_pheno_unrelated["phenoage_calibrated_advance"].std(axis=0)
df_pheno_unrelated["phenoage_calibrated_advance_std"] = (df_pheno_unrelated["phenoage_calibrated_advance"] - mean_val) / stdev_val

median_val = df_pheno_unrelated["phenoage_calibrated_advance"].median(axis=0)
lower_quartile = np.percentile(df_pheno_unrelated["phenoage_calibrated_advance"], 25)
upper_quartile = np.percentile(df_pheno_unrelated["phenoage_calibrated_advance"], 75)
min_val = df_pheno_unrelated["phenoage_calibrated_advance"].min(axis=0)
max_val = df_pheno_unrelated["phenoage_calibrated_advance"].max(axis=0)

df_pheno_unrelated_meta_added = df_pheno_unrelated.merge(df_meta, how="inner", on="IID")
sex_counts = df_pheno_unrelated_meta_added["Sex"].value_counts()

xticklabels = [
    f"Female\n(N={sex_counts.get('F', 0):,})",
    f"Male\n(N={sex_counts.get('M', 0):,})"
]

sns.boxplot(
    data=df_pheno_unrelated_meta_added,
    x="Sex",
    y="phenoage_calibrated_advance",
    width=0.5,
    color="grey",
    linewidth=0.8,
    ax=ax_3,
)

pairs = [("F", "M")]

annotator = Annotator(
    ax=ax_3,
    data=df_pheno_unrelated_meta_added,
    x="Sex",
    y="phenoage_calibrated_advance",
    order=["F", "M"],
    pairs=pairs
)

annotator.configure(test="Mann-Whitney", text_format="star", loc="inside", verbose=False)

annotator.apply_and_annotate()

ax_3.grid(axis="y", linewidth=0.5, linestyle="-")
ax_3.set_yticks(np.arange(-10, 90, 10))
ax_3.set_ylabel("PAA (calibrated)", fontsize=plt.rcParams["font.size"] + 6)
ax_3.set_xlabel("Sex", fontsize=plt.rcParams["font.size"] + 6)
ax_3.set_title(
    f"Korea10K\n(N={len(df_pheno_unrelated):,})",
    fontsize=plt.rcParams["font.size"] + 6
)

ax_3.set_xticklabels(xticklabels)

sns.despine(ax=ax_3)

path_assoc = "/BiO/Access/kyungwhan1998/genome/gwas/Results/KGP_creat_age_glucose_sbp_phenoage_calibrated_GWAS.phenoage_calibrated_advance_INT.glm.linear"
path_vep = "/BiO/Research/Korea10KGenome/Results/VariantAnnotation.by.VEP115/Korea10K_VariantAnnotation.FilterIndividualABHet2STD.VariantABHet0.4ABHom0.9ADSupp0.9.Biallelic.Geno_0.01.Mind_0.1.HWE_1e6.Het_3STD.Kinship_3rd.nonKorean.20251114.maf"
path_vep_gwas = "/BiO/Access/kyungwhan1998/genome/gwas/Results/KGP_creat_age_glucose_sbp_phenoage_calibrated_GWAS.phenoage_calibrated_advance_INT.glm.linear.annot"
gwas = pd.read_csv(path_assoc, sep="\t")
gwas_p_sorted = gwas[gwas["P"] < 1e-06].sort_values(by=["P"])

gwas_p_sorted["VAR"] = gwas_p_sorted["ID"] + "_" + gwas_p_sorted["REF"] + "_" + gwas_p_sorted["ALT"]
gwas_p_sorted["VAR_New"] = gwas_p_sorted["VAR"].apply(normalize_variant)

list_col_select = ['Chromosome',
 'vcf_pos',
 'Reference_Allele',
 'Allele1',
 'Allele2',
 'dbSNP_RS',
 'Hugo_Symbol']

with open(path_vep, mode="r") as fr:
    dict_col_select_ind = dict()
    header = fr.readline().rstrip("\n").split("\t")
    for ind, col in enumerate(header):
        if col in list_col_select:
            dict_col_select_ind[col] = ind
    
dict_col_select_ind_ref_allele1 = {k: v for k, v in dict_col_select_ind.items() if k != "Allele2"}
dict_col_select_ind_ref_allele2 = {k: v for k, v in dict_col_select_ind.items() if k != "Allele1"}

new_header = ['Chromosome',
 'vcf_pos',
 'Reference_Allele',
 'Alternate_Allele',
 'dbSNP_RS',
 'Hugo_Symbol']

if not os.path.exists(path_vep_gwas):
    with open(path_vep, mode="r") as fr , open(path_vep_gwas, mode="w") as fw:
        header = fr.readline().rstrip("\n").split("\t")
        write_header = "\t".join(new_header) + "\n"
        fw.write(write_header)
        
        for line in fr:
            record = line.rstrip("\n").split("\t")
            if record[dict_col_select_ind["Reference_Allele"]] == record[dict_col_select_ind["Allele1"]]:
                record_select_allele = list(map(lambda x: record[x], list(dict_col_select_ind_ref_allele2.values())))
                var_id = "_".join(record_select_allele[:4])
                if var_id in set(gwas_p_sorted["VAR_New"]):
                    write_content = "\t".join(record_select_allele) + "\n"
                    fw.write(write_content)
                    
            else:
                record_select_allele = list(map(lambda x: record[x], list(dict_col_select_ind_ref_allele1.values())))
                var_id = "_".join(record_select_allele[:4])
                if var_id in set(gwas_p_sorted["VAR_New"]):
                    write_content = "\t".join(record_select_allele) + "\n"
                    fw.write(write_content)

vep_gwas = pd.read_csv(path_vep_gwas, sep="\t")
vep_gwas["Annotation"] = vep_gwas["Hugo_Symbol"] + "\n" + "(" + vep_gwas["dbSNP_RS"] + ")"

ax, list_annot_text = manhattan_plot(gwas, sigline=1e-06, vep_df=vep_gwas, ax=ax_4)
ax_4.set_xlabel("Chromosome", fontsize=plt.rcParams["font.size"]+6)
ax_4.set_ylabel("$-log_{10}$(P)", fontsize=plt.rcParams["font.size"]+6)
ax_4.set_xticklabels(ax_4.get_xticklabels(), fontsize=plt.rcParams["font.size"]+1)
ax_4.set_yticklabels(ax_4.get_yticklabels(), fontsize=plt.rcParams["font.size"]+1)
ax_4.set_xmargin(0)
sns.despine(ax=ax_4, top=True, right=True)

lambda_median=calc_lambda_median_from_df(gwas) 
results = stats.probplot(gwas["P"], dist=stats.uniform, plot=None) 
ax_5.scatter(-np.log10(results[0][0]), -np.log10(results[0][1]), s=10, color="k") 
ax_5.plot([0, 8], [0, 8], color="firebrick") 
lambda_label = rf"$\lambda_{{GC}}$ = {float(lambda_median):.3f}" 
ax_5.text( 0.02, 0.98, 
          lambda_label, 
          transform=ax_5.transAxes, 
          fontsize=plt.rcParams["font.size"]+5, 
          va="top", ha="left", 
          bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8, edgecolor="none") ) 
ax_5.set_xlabel("Expected $-log_{10}$(P)", fontsize=plt.rcParams["font.size"]+5) 
ax_5.set_ylabel("Observed $-log_{10}$(P)", fontsize=plt.rcParams["font.size"]+5) 
sns.despine(ax=ax_5, top=True, right=True)

plt.show()
plt.close()

# %%
