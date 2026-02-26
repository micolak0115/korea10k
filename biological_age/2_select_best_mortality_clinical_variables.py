# %%
import glob
import os
import pickle

import pandas as pd

# %%
dict_clin_c_index = dict()
dir_cox_model = "/BiO/Research/GeroExpressome/Results/Cox_PH_Model"
clinical = "creat_waist_bmi_trig_age_glucose_sbp_hdl_lymph_dbp_totchol_monopa"
list_file_cox_models = glob.glob(f"{dir_cox_model}/cox_model*")
list_file_cox_models = list(filter(lambda x: clinical in x, list_file_cox_models))
list_clin = list(map(lambda x: os.path.basename(x).replace("cox_model_", "").replace(".pk", ""), list_file_cox_models))
for file_cox_model, clin in zip(list_file_cox_models, list_clin):
    with open(file_cox_model, mode="rb") as fr:
        cox_model = pickle.load(fr)
        dict_cox_params = dict(cox_model.params_)
        list_cox_params = list(dict_cox_params)
        num_cox_params = len(list_cox_params)
        dict_sig_hazard = dict(sorted({k: v for k, v in dict(cox_model.hazard_ratios_).items() if round(float(v), 2) != 1}.items(), key=lambda x: x[1])[::-1])
        list_sig_hazard = list(dict_sig_hazard.keys())
        list_hazard_ratio = list(dict_sig_hazard.values())
        num_sig_hazard = len(list_sig_hazard)
        c_index = round(float(cox_model.concordance_index_), 5)
        dict_clin_c_index[clin] = [list_sig_hazard, list_hazard_ratio, num_sig_hazard, c_index]

df_clin_c_index = pd.DataFrame.from_dict(dict_clin_c_index, orient="index")
df_clin_c_index.columns = ["Sig_Hazard", "Hazard_Ratio", "Num_Hazard", "C_index"]

df_clin_c_index_sort = df_clin_c_index.sort_values(by=["C_index"], ascending=False)

# df_clin_c_index_sort_filter_out_alb = df_clin_c_index_sort[df_clin_c_index_sort["Sig_Hazard"].apply(lambda x: "albumin" not in x)]
# df_clin_c_index_sort_filter_out_bmi = df_clin_c_index_sort_filter_out_alb[df_clin_c_index_sort_filter_out_alb["Sig_Hazard"].apply(lambda x: "bmi" not in x)]
# df_clin_c_index_sort_filter_out_creat = df_clin_c_index_sort[df_clin_c_index_sort["Sig_Hazard"].apply(lambda x: "creat" not in x)]

# %%
import matplotlib.pyplot as plt
import seaborn as sns

plt.figure(figsize=(1, 2))
ax = sns.heatmap(df_clin_c_index_sort[["C_index"]], cmap="rocket", vmin=0.88, vmax=0.90)
yticklabels = ax.get_yticklabels()
ax.set_yticklabels(yticklabels, fontsize=10, ha="right")

plt.show()
plt.close()

# %%
dict_hazard =  dict(cox_model.hazard_ratios_)
dict_univariate = {k: v for k, v in dict_hazard.items() if abs(round(v, 4)) > 1}
list_univariate = list(dict_univariate.keys())
# %%
