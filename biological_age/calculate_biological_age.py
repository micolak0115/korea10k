# %%
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# %%
path_phenoage = "/BiO/Research/GeroExpressome/Resources/Scripts/BiologicalRNAClock/KGP_creat_rdw_ttbl_hba1c_uap_monopa_bun_sbp_waist_wbc_dbp_lymph_rbc_age_phenoage_added.txt"
path_kgp_data = "/BiO/Research/GeroExpressome/Resources/Data/KOGIC/EHR/KU10K_CI_Korea10K_Raw_Ver2.7_incl_infectomics.tsv"

# %%
from matplotlib import font_manager as fm

path_font = "/BiO/Access/kyungwhan1998/miniconda3/lib/python3.11/site-packages/matplotlib/mpl-data/fonts/ttf/NanumGothic.ttf"
prop = fm.FontProperties(fname=path_font)
plt.rcParams['font.family'] = prop.get_name()

# %%
df_phenoage = pd.read_csv(path_phenoage, sep=",", index_col=["Sample_ID_New"])
df_kgp_data = pd.read_csv(path_kgp_data, sep="\t")

# %%
df_phenoage["age_acc"] = abs(df_phenoage["phenoage"] - df_phenoage["age"])
df_kgp_data["Sample_ID_New"] = df_kgp_data["Sample_ID"].astype(str) + "-" + df_kgp_data["Measurement_Number"].astype(str)
df_kgp_data["Sex"] = df_kgp_data["Sex"].apply(lambda x: 1 if x == "M" else 2)
kgp_data = df_kgp_data.set_index("Sample_ID_New")
df_phenoage["AMI"] = kgp_data["AMI"]
df_phenoage["Group"] = kgp_data["Group"]

import numpy as np
# %%
from scipy.stats import pearsonr
from sklearn.metrics import mean_absolute_error, r2_score

df_phenoage = df_phenoage.replace([np.inf, -np.inf], np.nan).dropna(subset=["age", "phenoage"])
r, p_value = pearsonr(df_phenoage["age"], df_phenoage["phenoage"])
r2 = r2_score(df_phenoage["age"], df_phenoage["phenoage"])
mae = mean_absolute_error(df_phenoage["age"], df_phenoage["phenoage"])

plt.figure(figsize=(5, 5))
sns.scatterplot(data=df_phenoage, x="age", y="phenoage", edgecolor="k", alpha=0.5, hue="Group", palette="Spectral")
# for idx, row in df_phenoage[df_phenoage["age_acc"] > 20].iterrows():
#     plt.text(row["age"], row["phenoage"], str(idx), fontsize=9, ha="left", va="bottom", color="k", weight="bold")
plt.plot([20, 90], [20, 90], color="k", linestyle="--")
# Add text with metrics
plt.text(70, 10, f"r = {r:.3f}\nR² = {r2:.3f}\nMAE = {mae:.3f}", 
         fontsize=12, weight="bold")
plt.legend("", frameon=False)
plt.show()
plt.close()

# %%
path_font = "/BiO/Access/kyungwhan1998/miniconda3/lib/python3.11/site-packages/matplotlib/mpl-data/fonts/ttf/Arial.ttf"
prop = fm.FontProperties(fname=path_font)
plt.rcParams['font.family'] = prop.get_name()

path_phenoage = "/BiO/Research/GeroExpressome/Resources/Scripts/BiologicalRNAClock/KGP_creat_rdw_ttbl_hba1c_uap_monopa_bun_sbp_waist_wbc_dbp_lymph_rbc_age_phenoage_added.txt"
path_ehr_inf = "/BiO/Research/GeroExpressome/Resources/Data/KOGIC/EHR/infectomics_CRF_20230410_edit.xlsx"

df_phenoage = pd.read_csv(path_phenoage, sep=",")
df_phenoage["Sample_ID"] = df_phenoage["Sample_ID_New"].apply(lambda x: "-".join(x.split("-")[:-1]))
df_ehr_inf = pd.read_excel(path_ehr_inf, skiprows=1)
df_severity = df_ehr_inf[["Subject NO.(고유번호)", "중증도분류"]]
df_severity.columns = ["Sample_ID", "Severity"]
df_merged = df_phenoage.merge(df_severity, how="inner", on="Sample_ID")

plt.figure(figsize=(6, 6))
sns.scatterplot(df_merged, x="Severity", y="phenoage_advance", color="grey")
for idx, row in df_merged[df_merged["phenoage_advance"] > 5].iterrows():
    dict_row = dict(row)
    sample_id = dict_row["Sample_ID_New"]
    plt.text(row["Severity"], row["phenoage_advance"], str(sample_id), fontsize=9, ha="left", va="bottom", color="k", weight="bold")
plt.xticks(range(1, len(df_merged["Severity"].unique())+1))
plt.xlabel("Severity", fontsize=16)
plt.ylabel("Age Acceleration", fontsize=16)
plt.grid(axis="both")


# %%
import seaborn as sns
from matplotlib.patches import Rectangle
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
from sklearn.linear_model import BayesianRidge

# %%
kgp_data_clean = pd.read_csv("/BiO/Research/GeroExpressome/Resources/Scripts/BiologicalRNAClock/KGP_biomarkers_only_common_NHANES3.txt", sep="\t")
kgp_data_clean_no_sample_id = kgp_data_clean.drop(columns=["Sample_ID_New"])
kgp_data_bm_only = kgp_data_clean_no_sample_id[sig_biomarkers]
bayesian_ridge_estimator = BayesianRidge()
mice_imputer = IterativeImputer(estimator=bayesian_ridge_estimator, max_iter=10, random_state=0)
kgp_data_imputed = pd.DataFrame(mice_imputer.fit_transform(kgp_data_bm_only), columns=kgp_data_bm_only.columns)

# %%
dict_marker_num_na = dict()
for i, bm in enumerate(biomarkers):
    if bm in kgp_data_clean.columns:
        marker_na = list(filter(lambda x: x=="nan", kgp_data_clean[bm].astype(str)))
        num_na = len(marker_na)
        dict_marker_num_na[bm] = num_na
dict_marker_num_na_sort = dict(sorted(dict_marker_num_na.items(), key=lambda x:x[1]))
plt.bar(x=dict_marker_num_na_sort.keys(), height=dict_marker_num_na_sort.values(), color="grey")
plt.xticks(fontsize=12, rotation=90)
plt.xlabel("Clinical Values", fontsize=15)
plt.ylabel("Number of Missing Data", fontsize=15)

# %%
bayesian_ridge_estimator = BayesianRidge()
mice_imputer = IterativeImputer(estimator=bayesian_ridge_estimator, max_iter=10, random_state=0)
kgp_data_imputed = pd.DataFrame(mice_imputer.fit_transform(kgp_data), columns=kgp_data.columns)

kgp_data["source"] = "original"
kgp_data_imputed["source"] = "imputed"

data_merged = pd.concat([kgp_data, kgp_data_imputed], axis=0).reset_index(drop=True)
num_columns = data_merged.select_dtypes(include=["number"]).columns.tolist()
num_plots = len(num_columns)
fig, axes = plt.subplots(nrows=(num_plots // 2) + (num_plots % 2), ncols=5, figsize=(12, num_plots*1.5))
axes = axes.flatten()

custom_palette = {"imputed": "firebrick", "original": "grey"}
for i, col in enumerate(num_columns):
    ax = axes[i]
    sns.boxplot(x="source", y=col, data=data_merged, ax=ax, hue="source", order=list(custom_palette.keys()), palette=custom_palette)
    ax.set_title(col)
    ax.set_xlabel("")
    ax.set_xticklabels(labels=list(custom_palette.keys()), rotation=45, rotation_mode="anchor", ha="right")

for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])

handles = [Rectangle((0, 0), 1, 1, color=c, ec="k") for c in list(custom_palette.values())]
fig.legend(handles=handles, labels=list(custom_palette.keys()), loc="center left", bbox_to_anchor=(0.9, 0.75), fontsize=14, frameon=False)

plt.tight_layout(rect=[0, 0, 0.9, 1])
plt.show()
plt.close()
