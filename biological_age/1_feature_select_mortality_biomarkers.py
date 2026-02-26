# %%
import os
import subprocess

import pandas as pd

import mortality
from job_maker import write_script

# %%
INDIR = "/BiO/Research/GeroExpressome/Resources/Data/NHANES"
WORKDIR = "/BiO/Research/GeroExpressome/Resources/Scripts/MortalityClock"
OUTDIR = "/BiO/Research/GeroExpressome/Results"
os.makedirs(os.path.join(OUTDIR, "Cox_PH_Model"), exist_ok=True)
nhanes3_data = pd.read_csv(f"{INDIR}/NHANES3.txt", sep=",")
nhanes4_data = pd.read_csv(f"{INDIR}/NHANES4.txt", sep=",")
var_conv_json = f"{WORKDIR}/KU10K_NHANES_clinical_variable_conversion_table.json"
path_rm_biomarker = f"{WORKDIR}/remove_biomarkers_final.txt"
path_nhanes3_data_biomarkers_only = f"{OUTDIR}/NHANES3_biomarkers_only_final.txt"
dir_script = os.path.join(OUTDIR, "Script")
option = "min_plus_std"
n_threads = 1
n_jobs = 7

# %%
nhanes3_data_clean = mortality.parse_nhanes_data_mortality_stats(nhanes3_data)
nhanes3_data_clean["gender"] = nhanes3_data_clean["gender"].apply(lambda x: 0 if x == 1 else 1)
remove_biomarkers = mortality.get_remove_biomakers(path_rm_biomarker)
var_conv_table = mortality.get_conversion_table_clinical_variables(var_conv_json)
usable_biomarkers = mortality.get_selected_biomarkers(var_conv_table, remove_biomarkers, list_manual=[])
mortality.get_nhanes_selected_biomarker(nhanes3_data_clean, usable_biomarkers, outfile=path_nhanes3_data_biomarkers_only)

# %%
@write_script
def job(data, biomarkers, workdir, outdir, option, n_threads, n_jobs, job_name):
    job_cmd = f"python {workdir}/cox_fitter.py --data {data} --biomarkers {biomarkers} --outdir {outdir} --option {option} --threads {n_threads}"
    print(f"running {job_name} with {n_jobs} jobs")
    
    return job_cmd

job_name = "job_name_" + "_".join(usable_biomarkers)
job(data=path_nhanes3_data_biomarkers_only, biomarkers=" ".join(usable_biomarkers), workdir=WORKDIR, outdir=OUTDIR, option=option, n_threads=n_threads, n_jobs=n_jobs, job_name=job_name)
path_script = os.path.join(dir_script, f"{job_name}.sh")
subprocess.run(f"qsub {path_script}", shell=True)

# %%
