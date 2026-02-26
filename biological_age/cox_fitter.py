import argparse
import os
import pickle

import pandas as pd

import mortality


# %%
def run(path_data, biomarkers, outdir, option, threads):
    data = pd.read_csv(path_data, sep="\t")
    biomarkers_for_filename = "_".join(biomarkers)
    path_dict_params_cox_model = f"{outdir}/Cox_PH_Model/dict_params_{biomarkers_for_filename}.pk"
    data_bm = mortality.get_nhanes_selected_biomarker(data, biomarkers, outfile=None)
    dict_param = mortality.run_kfold_cv_select_reg_param_cox_PH_model(data, biomarkers, path_dict_params_cox_model, n_jobs=threads)
    mortality.plot_partial_LLH_deviance(dict_param, biomarkers, outdir)
    penalizer = mortality.get_penalizer_optimum(dict_param, option=option)
    cox_model = mortality.train_cox_PH_model(data_bm, penalizer)
    outmodel = f"{outdir}/Cox_PH_Model/cox_model_{biomarkers_for_filename}_{option}_deviance.pk"
    with open(outmodel, 'wb') as handle:
        pickle.dump(cox_model, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--data")
    parser.add_argument("-b", "--biomarkers", nargs="+")
    parser.add_argument("-od", "--outdir")
    parser.add_argument("-op", "--option", default="min")
    parser.add_argument("-t", "--threads", default=1, type=int)
    args = parser.parse_args()
    dict_args = vars(args)
    os.makedirs(dict_args["outdir"], exist_ok=True) 
    run(dict_args["data"], dict_args["biomarkers"], dict_args["outdir"], option=dict_args["option"], threads=dict_args["threads"])