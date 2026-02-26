
# %%
import json
import math
import os
import pickle

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from lifelines import CoxPHFitter
from scipy.interpolate import interp1d
from scipy.stats import ranksums
from statannotations.Annotator import Annotator


# %%
def parse_nhanes_data_mortality_stats(nhanes_data, ucod_leading_incl=[4, 8, 10]):
    nhanes_data_eligstat_filt = nhanes_data[nhanes_data["eligstat"] == 1]
    cond_nan_ucod_leading_but_dead = np.logical_and(nhanes_data_eligstat_filt["ucod_leading"].isna(), nhanes_data_eligstat_filt["status"]==1)
    cond_ucod_leading_but_alive = np.logical_and(nhanes_data_eligstat_filt["ucod_leading"].notna(), nhanes_data_eligstat_filt["status"]==0)
    cond_ucod_issued = np.logical_or(cond_nan_ucod_leading_but_dead, cond_ucod_leading_but_alive)
    nhanes_data_eligstat_ucod_issue_filt = nhanes_data_eligstat_filt[~cond_ucod_issued]
    cnd_ucod_excl = np.logical_and(nhanes_data_eligstat_ucod_issue_filt["ucod_leading"].notna(), nhanes_data_eligstat_ucod_issue_filt["ucod_leading"].isin(ucod_leading_incl))
    nhanes_data_ucod_select = nhanes_data_eligstat_ucod_issue_filt[~cnd_ucod_excl]
    nhanes_data_clean = nhanes_data_ucod_select[nhanes_data_ucod_select["time"]!=0]

    return nhanes_data_clean

def plot_follow_up_period(data, step=60):
    follow_up_period = data[data["status"] == 1]["time"]
    plt.hist(follow_up_period, color="grey")
    for i in range(0, int(max(follow_up_period)+1), step):
        plt.axvline(i, color="firebrick", linestyle="--")
    plt.title("F/U Period for Deceased")
    plt.xlabel("Examination~Death(Months)")
    plt.ylabel("Frequency")
    
def get_conversion_table_clinical_variables(var_conv_json):
    with open(var_conv_json, mode="rb") as frb:
        var_conv_table = json.load(frb)
    
    return var_conv_table

def get_remove_biomakers(path_rm_biomarkers):
    with open(path_rm_biomarkers, mode="r") as fr:
        remove_biomarkers = fr.readlines()
    remove_biomarkers = list(map(lambda x: x.rstrip("\n"), remove_biomarkers))
    
    return remove_biomarkers

def get_selected_biomarkers(var_conv_table, remove_biomarkers, list_manual=["lncrp"]):
    list_cols = list(var_conv_table.keys())
    biomarkers = list(set(list_cols).difference(set(remove_biomarkers)))
    for manual in list_manual:
        if manual not in biomarkers:
            biomarkers.append(manual)
    
    return biomarkers

def plot_hist_num_na_biomarkers(data, biomarkers):
    dict_marker_num_na = dict()
    for i, bm in enumerate(biomarkers):
        if bm in data.columns:
            marker_na = list(filter(lambda x: x=="nan", data[bm].astype(str)))
            num_na = len(marker_na)
            dict_marker_num_na[bm] = num_na
    dict_marker_num_na_sort = dict(sorted(dict_marker_num_na.items(), key=lambda x:x[1]))
    plt.bar(x=dict_marker_num_na_sort.keys(), height=dict_marker_num_na_sort.values(), color="grey")
    plt.xticks(fontsize=12, rotation=90)
    plt.xlabel("Clinical Values", fontsize=15)
    plt.ylabel("Number of Missing Data", fontsize=15)

def get_nhanes_selected_biomarker(data, biomarkers, outfile=None):
    data_bm = data[biomarkers + ["time", "status"]].dropna()
    if outfile is None:
        pass
 
    else:
        data_bm.to_csv(outfile, sep="\t", index=False)
    
    return data_bm

def run_kfold_cv_select_reg_param_cox_PH_model(data, biomarkers, outfile, col_duration = "time", col_event = "status", l1_ratio=0.5, cv=10, random_state=42, n_jobs=10):
    from run_lifelines_cv import select_best_penalizer
    data_bm = get_nhanes_selected_biomarker(data, biomarkers)
    table = data_bm.reset_index(drop=True)
    penalizer_range = 10 ** np.linspace(-4, 1, 41)
    dict_param = select_best_penalizer(table, col_duration, col_event, penalizer_range, l1_ratio, cv = cv, random_state = random_state, n_jobs = n_jobs)
    with open(outfile, 'wb') as handle:
        pickle.dump(dict_param, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    return dict_param

def plot_partial_LLH_deviance(dict_param, biomarkers, outdir):
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = 14
    plt.figure(figsize = (9, 5))
    plt.boxplot(dict_param["Partial_LLH_deviance"], patch_artist = True, boxprops=dict(facecolor='w'))
    plt.xticks(range(1, len(dict_param["Penalizers"])+1), list(map(lambda val: f"{val:.2f}", np.log10(dict_param["Penalizers"]))), zorder = 3, fontsize = plt.rcParams["font.size"]-3, rotation = 90)
    plt.axvline(dict_param["Index_min_deviance"]+1, linestyle = "--", linewidth = 1, color = "gray", zorder = 1, label = f"Minimum Deviance\n(Penalizer = {dict_param['Penalizers'][dict_param['Index_min_deviance']]:.3f})")
    plt.axvline(dict_param["Index_min_plus_std_deviance"]+1, linestyle = "--", linewidth = 1, color = "cornflowerblue", zorder = 1, label = f"Minimum+1STD Deviance\n(Penalizer = {dict_param['Penalizers'][dict_param['Index_min_plus_std_deviance']]:.3f})")
    plt.xlabel("Log10(Penalizer)")
    plt.ylabel("Partial Log Likelihood Deviance")
    plt.legend(loc = "upper left", bbox_to_anchor = (0.01, 0.99), fontsize = plt.rcParams["font.size"]-1)
    plt.tight_layout()
    biomarkers_for_filename = "_".join(biomarkers)
    outfile = os.path.join(outdir, "Cox_PH_Model", f"plot_partial_LLH_deviance_{biomarkers_for_filename}.png")
    plt.savefig(outfile, dpi=300)

def get_penalizer_optimum(dict_param, option):
    if option == "min":
        penalizer = dict_param['Penalizers'][dict_param["Index_min_deviance"]]

    elif option == "min_plus_std":
        penalizer = dict_param['Penalizers'][dict_param["Index_min_plus_std_deviance"]]

    elif option == "min_plus_half_std":
        penalizer = dict_param['Penalizers'][dict_param["Index_min_plus_half_std_deviance"]]
    
    else:
        raise Exception
    
    return penalizer

def train_cox_PH_model(data_bm, penalizer, l1_ratio=0.5, duration_col='time', event_col='status', verbose=True):
    cox_model = CoxPHFitter(penalizer=penalizer, l1_ratio=l1_ratio)
    cox_model.fit(data_bm, duration_col=duration_col, event_col=event_col)
    if verbose:
        cox_model.print_summary()
    
    return cox_model

def get_biomarkers_hazard(cox_model):
    dict_hazard = cox_model.hazard_ratios_.to_dict()
    dict_hazard_sig = {k: v for k, v in dict_hazard.items() if round(float(v), 5)}
    
    return dict_hazard_sig

def attach_baseline_and_predicted_hazard(cox_model, df_partial, df_source, col_age="age"):
    # --- Interpolate baseline hazard ---
    bh = cox_model.baseline_hazard_.reset_index()
    bh.columns = ["time", "baseline_hazard"]
    interp_baseline = interp1d(
        bh["time"].values,
        bh["baseline_hazard"].values,
        kind="linear",
        bounds_error=False,
        fill_value="extrapolate"
    )
    
    common_ids = df_partial.index.intersection(df_source.index)

    df = df_partial.loc[common_ids].copy()
    df["age"] = df_source.loc[common_ids, col_age].values
    # Interpolate baseline hazard for each sample's age
    df["baseline_hazard"] = interp_baseline(df["age"])
    # Convert log10(partial_hazard) to actual partial hazard
    df["partial_hazard"] = 10 ** df["log10(partial_hazard)"]
    # Compute predicted absolute hazard
    df["predicted_hazard"] = df["baseline_hazard"] * df["partial_hazard"]
    # Difference from baseline
    df["hazard_over_baseline"] = df["predicted_hazard"] / df["baseline_hazard"]
    df["log10(hazard_over_baseline)"] = df["hazard_over_baseline"].apply(lambda x: math.log10(x))
    
    return df

def calculate_partial_hazard(cox_model, data, biomarkers, col_status=None, col_age = "age", col_sample="Sample_ID_New", col_ph="partial_hazard"):
    if col_status is None:
        test_data_bm_only = data[biomarkers + [col_sample]].dropna()
    else:
        test_data_bm_only = data[biomarkers + [col_status, col_sample]].dropna()
    test_ph = cox_model.predict_partial_hazard(test_data_bm_only)
    test_data_bm_only_predictions_attached = pd.concat([test_data_bm_only, test_ph], axis=1)
    test_data_bm_only_predictions_attached = test_data_bm_only_predictions_attached.rename(columns={0: col_ph})
    test_data_bm_only_predictions_attached[f"log10({col_ph})"] = test_data_bm_only_predictions_attached[col_ph].apply(math.log10)
    test_data_bm_only_predictions_and_baselline_attached = attach_baseline_and_predicted_hazard(cox_model, test_data_bm_only_predictions_attached, data, col_age=col_age)
    
    return test_data_bm_only_predictions_and_baselline_attached

def plot_partial_hazard(test_data_bm_only_predictions_attached, event_col='status', col_ph="partial_hazard"):
    ax = sns.boxplot(
        data=test_data_bm_only_predictions_attached, 
        x=event_col, 
        y=col_ph, 
        color="grey",
        width=0.5
    )

    # --- Add annotation using statannotations ---
    pairs = [(0, 1)]
    annotator = Annotator(ax, pairs, data=test_data_bm_only_predictions_attached, x=event_col, y=col_ph)
    annotator.configure(test='Mann-Whitney', text_format='star', loc='inside', verbose=0)
    annotator.apply_and_annotate()

    # --- Labeling ---
    ax.set_ylabel("$log_{10}$(Partial Hazard)", fontsize=15)
    ax.set_xlabel("Death Status", fontsize=15)
    ax.set_xticks([0, 1], ["Alive", "Deceased"])
    
    return ax
# %%
