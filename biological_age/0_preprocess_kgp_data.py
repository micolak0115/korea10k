# %%
import json
import math
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

# %%
nhanes_data = pd.read_csv("/BiO/Research/GeroExpressome/Resources/Data/NHANES/NHANES3.txt", sep=",")

# %%
def remove_visit_number(sample_id):
    pattern = r"-V\d+"
    sample_id_rmv_visit = re.sub(pattern, "", sample_id)
    return sample_id_rmv_visit

# %%
kgp_data = pd.read_csv("/BiO/Research/GeroExpressome/Resources/Data/KOGIC/EHR/KU10K_CI_Korea10K_Raw_Ver2.7_incl_infectomics.tsv", sep="\t")
kgp_data["Sample_ID_New"] = kgp_data["Sample_ID"].astype(str) + "-" + kgp_data["Measurement_Number"].astype(str)
kgp_data["Sample_ID_New"] = kgp_data["Sample_ID_New"].apply(remove_visit_number)
kgp_data["Sex"] = kgp_data["Sex"].apply(lambda x: 1 if x == "M" else 2)
kgp_data = kgp_data.set_index("Sample_ID_New")

# %%
from collections import Counter

plt.rcParams["font.size"] = 15

measure_num = kgp_data[kgp_data["Age"].isnull()]["Measurement_Number"].to_list()
dict_measure_num = dict(Counter(measure_num))
plt.bar(dict_measure_num.keys(), dict_measure_num.values(), color="grey")
plt.xlabel("Measurement_Number")
plt.ylabel("Count")
plt.show()
plt.close()

# %%
kgp_data_age_notna = kgp_data[~kgp_data["Age"].isna()]
plt.hist(kgp_data_age_notna["Age"], color="grey")
plt.axvline(x=np.mean(kgp_data_age_notna["Age"]), color="k", linestyle="--")
plt.text(np.mean(kgp_data_age_notna["Age"])+1, 1500, s=f"mean age = {int(np.mean(kgp_data_age_notna['Age']))}", weight="bold")
plt.xlabel("Age")
plt.ylabel("Count")
plt.show()
plt.close()

# %%
dict_cnt_sex = dict(Counter(kgp_data_age_notna["Sex"]))
plt.bar(dict_cnt_sex.keys(), dict_cnt_sex.values(), color="grey")
plt.xticks(range(1, len(dict_cnt_sex.keys())+1))
plt.xlabel("Sex")
plt.ylabel("Count")
plt.show()
plt.close()

# %%
def remove_special_characters(x):
    cleaned_x = re.sub(r"[^\d.-]", "", str(x)).strip()
    if cleaned_x == "":
        return np.nan
    
    return float(cleaned_x)

# %%
dict_ind_val_crp = dict()
for ind, crp in kgp_data_age_notna[["hs-CRP", "CRP_Quantitative"]].iterrows():
    row = dict(crp)
    val_crps = list(row.values())
    val_crps_nona = list(filter(lambda x: str(x) != "nan", list(val_crps)))
    if len(val_crps_nona) == 0:
        val_crp = np.nan
    else:
        val_crp = val_crps_nona[0]
        val_crp = remove_special_characters(val_crp)
    dict_ind_val_crp[ind] = val_crp

kgp_data_age_notna["CRP"] = list(dict_ind_val_crp.values())
kgp_data_age_notna["lnCRP"] = kgp_data_age_notna["CRP"].apply(lambda x: math.log(x+0.001))

# %%
remove_biomarkers = ["cadmium", "insulin", "gender"]
var_conv_json = "/BiO/Research/GeroExpressome/Resources/Data/KOGIC/EHR/KU10K_NHANES_clinical_variable_conversion_table.json"
with open(var_conv_json, mode="rb") as frb:
    var_conv_table = json.load(frb)
var_conv_table_rev = {v: k for k,v in var_conv_table.items()}
list_cols = list(var_conv_table_rev.values())
biomarkers = list(set(list_cols).difference(set(remove_biomarkers)))
kgp_data_selected = kgp_data_age_notna[list(var_conv_table_rev.keys())]
kgp_data_selected.columns = list(var_conv_table_rev.values())
kgp_data_filtered = kgp_data_selected[biomarkers]
nhanes_data_filtered = nhanes_data[biomarkers]
kgp_data_cleaned = kgp_data_filtered.applymap(remove_special_characters)

# %%
nhanes_data_cleaned = nhanes_data_filtered.applymap(remove_special_characters)

# %%
dict_marker_num_na = dict()
for i, bm in enumerate(biomarkers):
    if bm in nhanes_data_cleaned.columns:
        marker_na = list(filter(lambda x: x=="nan", nhanes_data_cleaned[bm].astype(str)))
        num_na = len(marker_na)
        dict_marker_num_na[bm] = num_na
dict_marker_num_na_sort = dict(sorted(dict_marker_num_na.items(), key=lambda x:x[1]))

plt.figure(figsize=(12, 4))
plt.bar(x=dict_marker_num_na_sort.keys(), height=dict_marker_num_na_sort.values(), color="grey")
plt.xticks(rotation=45, rotation_mode="anchor", ha="right", fontsize=13)
plt.xlabel("Clinical Values", fontsize=15)
plt.ylabel("Number of Missing Data", fontsize=15)
plt.ylim(0, 20000)
plt.title("NHANES III", fontsize=16)
plt.show()
plt.close()

# %%
dict_marker_num_na = dict()
for i, bm in enumerate(biomarkers):
    if bm in kgp_data_cleaned.columns:
        marker_na = list(filter(lambda x: x=="nan", kgp_data_cleaned[bm].astype(str)))
        num_na = len(marker_na)
        dict_marker_num_na[bm] = num_na
dict_marker_num_na_sort = dict(sorted(dict_marker_num_na.items(), key=lambda x:x[1]))

plt.figure(figsize=(12, 4))
plt.bar(x=dict_marker_num_na_sort.keys(), height=dict_marker_num_na_sort.values(), color="grey")
plt.xticks(rotation=45, rotation_mode="anchor", ha="right", fontsize=13)
plt.xlabel("Clinical Values", fontsize=15)
plt.ylabel("Number of Missing Data", fontsize=15)
plt.ylim(0, 20000)
plt.title(f"Korean Genome Project", fontsize=16)
plt.show()
plt.close()

# %%
remove_biomarkers_based_on_missingness = ["fm", "ffm", "basopa", "eosnpa", "ldl", "ggt"]
biomarkers = list(set(kgp_data_cleaned.columns).difference(set(remove_biomarkers_based_on_missingness)))
kgp_data_cleaned_bm_only = kgp_data_cleaned[biomarkers]
path_kgp_bm_only = "/BiO/Research/GeroExpressome/Resources/Scripts/MortalityClock/KGP_biomarkers_only_common_NHANES3.txt"
kgp_data_cleaned_bm_only.to_csv(path_kgp_bm_only, sep="\t", index=True)

# %%
path_nhanes_bm_only = "/BiO/Research/GeroExpressome/Resources/Scripts/MortalityClock/NHANES3_biomarkers_only_common_KGP.txt"
biomarkers = list(set(kgp_data_cleaned.columns).difference(set(remove_biomarkers_based_on_missingness)))
nhanes_data_cleaned_bm_only = nhanes_data_cleaned[biomarkers]
nhanes_data_cleaned_bm_only.to_csv(path_nhanes_bm_only)

# %%
data = [nhanes_data_cleaned_bm_only, kgp_data_cleaned_bm_only]
colors = ["firebrick", "grey"]
labels = ["NHANES", "KGP"]

nrows = int(np.sqrt(len(biomarkers)))+1
ncols = int(np.sqrt(len(biomarkers)))+1
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(16, 10))
axes = axes.flatten()

# Plotting histograms
for i, bm in enumerate(biomarkers):
    for dat, color in zip(data, colors):
        if bm in dat.columns:
            marker_val = list(map(float, dat[bm].dropna()))
            mean_val = np.mean(marker_val)
            # std_val = np.std(marker_val)
            axes[i].hist(marker_val, bins=20, color=color, alpha=0.7)
            axes[i].set_title(bm, fontsize=16)
            axes[i].set_xlabel("Value", fontsize=14)
            axes[i].set_ylabel("Frequency", fontsize=14)
            handles = [Rectangle((0, 0), 1, 1, color=c, ec="k") for c in colors]
            axes[i].axvline(x=mean_val, color=color, linestyle="--")

        else:
            axes[i].axis("off")

for j in range(i + 1, len(axes)):
    axes[j].axis("off")

handles = [Rectangle((0, 0), 1, 1, color=c, ec="k") for c in colors]
fig.legend(handles=handles, labels=labels, loc="center left", bbox_to_anchor=(0.9, 0.5), fontsize=14, frameon=False)

plt.tight_layout(rect=[0, 0, 0.9, 1])

plt.show()
plt.close()
# %%
import seaborn as sns
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import fdrcorrection

nhanes_data_cleaned_bm_only["source"] = "NHANES"
kgp_data_cleaned_bm_only["source"] = "KGP"
data_merged = pd.concat([nhanes_data_cleaned_bm_only, kgp_data_cleaned_bm_only], axis=0).reset_index(drop=True)

num_columns = data_merged.select_dtypes(include=["number"]).columns.tolist()
p_values = []
for col in num_columns:
    group1 = data_merged[data_merged["source"] == "NHANES"][col].dropna()
    group1_rmv_chr = list(map(lambda x: remove_special_characters(x), group1))
    group2 = data_merged[data_merged["source"] == "KGP"][col].dropna()
    group2_rmv_chr = list(map(lambda x: remove_special_characters(x), group2))
    _, p_value = ttest_ind(group1_rmv_chr, group2_rmv_chr, equal_var=False)
    p_values.append(p_value)

_, fdr_values = fdrcorrection(p_values, alpha=0.05)


num_plots = len(num_columns)
fig, axes = plt.subplots(nrows=(num_plots // 2) + (num_plots % 2), ncols=10, figsize=(20, num_plots*2))
axes = axes.flatten()

custom_palette = {"NHANES": "firebrick", "KGP": "grey"}
for i, (col, fdr) in enumerate(zip(num_columns, fdr_values)):
    ax = axes[i]
    sns.boxplot(x="source", y=col, data=data_merged, ax=ax, hue="source", palette=custom_palette)
    fdr_text = f"FDR = {fdr:.1E}"
    color = "red" if fdr < 0.05 else "black"
    ax.set_title(col, color=color)
    ax.set_xticklabels(labels=list(custom_palette.keys()), rotation=90)
    ax.set_xlabel("")

for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])

handles = [Rectangle((0, 0), 1, 1, color=c, ec="k") for c in colors]
fig.legend(handles=handles, labels=labels, loc="center left", bbox_to_anchor=(0.9, 0.75), frameon=False)

plt.subplots_adjust(wspace=0.2, hspace=0.2)
plt.tight_layout(rect=[0, 0, 0.9, 1])
plt.show()
plt.close()

# %%
list_ku10k = list(filter(lambda x: "U10K" in x, kgp_data_cleaned_bm_only.index))
kgp_data_cleaned_bm_only_ku10k = kgp_data_cleaned_bm_only.loc[list_ku10k]

dict_bm_num_avail = dict()
for bm in kgp_data_cleaned_bm_only_ku10k.columns:
    num_all = len(kgp_data_cleaned_bm_only_ku10k)
    num_na = sum(kgp_data_cleaned_bm_only_ku10k[bm].isna())
    num_avail = num_all - num_na
    dict_bm_num_avail[bm] = num_avail

dict_bm_num_avail_sort = dict(sorted(dict_bm_num_avail.items(), key=lambda x: -x[1]))

list_rm_biomarkers = list()
for bm, cnt in dict_bm_num_avail_sort.items():
    if cnt < 5000:
        list_rm_biomarkers.append(bm)

list_rm_biomarkers

# %%
plt.figure(figsize=(10, 5))
plt.bar(dict_bm_num_avail_sort.keys(), dict_bm_num_avail_sort.values(), color="grey")
plt.xlabel("Biomarkers", fontsize=20)
plt.ylabel("Sample Count", fontsize=20)
plt.xticks(list(dict_bm_num_avail_sort.keys()), rotation=45, rotation_mode="anchor", ha="right")
plt.axhline(y=5000, linestyle="--", color="firebrick", linewidth=2)
plt.margins(x=0)
plt.show()
plt.close()

# %%
list_covid = list(filter(lambda x: str(x).startswith("C19"), list(kgp_data_cleaned_bm_only.index)))
dict_bm_num_avail = dict()
for bm in kgp_data_cleaned_bm_only.loc[list_covid]:
    num_all = len(kgp_data_cleaned_bm_only.loc[list_covid])
    num_na = sum(kgp_data_cleaned_bm_only.loc[list_covid][bm].isna())
    num_avail = num_all - num_na
    dict_bm_num_avail[bm] = num_avail

rmv_biomarkers_for_c19 = list()
for bm, num_miss in dict_bm_num_avail.items():
    if int(num_miss) < 150:
        rmv_biomarkers_for_c19.append(bm)

rmv_biomarkers_for_c19
# %%
