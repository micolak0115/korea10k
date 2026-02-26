# %%
import math
import os
import pickle

import numpy as np
import pandas as pd

import mortality

# %%
option = "min_plus_std"
# clinical = "creat_waist_bmi_trig_age_glucose_sbp_hdl_lymph_dbp_totchol_monopa"
clinical = "age_wbc_ttbl_alp_monopa_lymph_albumin_glucose_uap_bun_rbc_creat"

path_cox_model = f"/BiO/Research/GeroExpressome/Results/Cox_PH_Model/cox_model_{clinical}_{option}_deviance.pk"
outdir = "/BiO/Research/GeroExpressome/Results/Cox_PH_Model"

path_nhanes_iii = "/BiO/Research/GeroExpressome/Resources/Data/NHANES/NHANES3.txt"
path_nhanes_iv = "/BiO/Research/GeroExpressome/Resources/Data/NHANES/NHANES4.txt"
path_kgp = f"/BiO/Research/GeroExpressome/Results/PreprocessedData/KGP_biomarkers_only_common_NHANES3.txt"

outfilename = "log_transformed_partial_hazard_{0}_{1}.tsv"

# %%
with open(path_cox_model, mode="rb") as fb:
    cox_model = pickle.load(fb)

# %%
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

# %%
dict_biomarkers_hazard = mortality.get_biomarkers_hazard(cox_model)
biomarkers = list(dict_biomarkers_hazard.keys())

df_partial_hazard_nhanes_iii = mortality.calculate_partial_hazard(
    cox_model, nhanes_iii_processed_data.reset_index(), biomarkers, col_status="status"
).set_index("Sample_ID_New")
df_partial_hazard_nhanes_iii.to_csv(os.path.join(outdir, outfilename.format("NHANESIII", clinical)),
                                   sep="\t", index=True)

df_partial_hazard_nhanes_iv = mortality.calculate_partial_hazard(
    cox_model, nhanes_iv_processed_data.reset_index(), biomarkers, col_status="status"
).set_index("Sample_ID_New")
df_partial_hazard_nhanes_iv.to_csv(os.path.join(outdir, outfilename.format("NHANESIV", clinical)),
                                   sep="\t", index=True)

df_partial_hazard_kgp = mortality.calculate_partial_hazard(
    cox_model, kgp_data.reset_index(), biomarkers
).set_index("Sample_ID_New")
df_partial_hazard_kgp.to_csv(os.path.join(outdir, outfilename.format("KGP", clinical)),
                             sep="\t", index=True)

# %%
ax_0 = mortality.plot_partial_hazard(df_partial_hazard_nhanes_iii, col_ph="log10(partial_hazard)")
ax_0.set_title(f"NHANESIII (N={len(df_partial_hazard_nhanes_iii)})")
ax_1 = mortality.plot_partial_hazard(df_partial_hazard_nhanes_iv, col_ph="log10(partial_hazard)")
ax_1.set_title(f"NHANESIV (N={len(df_partial_hazard_nhanes_iv)})")
