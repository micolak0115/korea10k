# %%
import os
import subprocess


# %%
def get_cmd(bfile, covar, covar_name, outprefix, pheno, pheno_name):
    cmd = f'''#!bin/sh
    #$ -N job.Run_GWAS
    #$ -o /BiO/Access/kyungwhan1998/genome/gwas/Results/QSUB/job.Run_GWAS.out
    #$ -e /BiO/Access/kyungwhan1998/genome/gwas/Results/QSUB/job.Run_GWAS.err
    #$ -S /bin/bash
    #$ -cwd
    #$ -m bea
    #$ -pe smp 2
    echo "=========================================================="
    echo "Starting on       : $(date)"
    echo "Running on node   : $(hostname)"
    echo "Current directory : $(pwd)"
    echo "Current job ID    : $JOB_ID"
    echo "Current job name  : $JOB_NAME"
    #echo "Task index number : $TASK_ID"
    echo "=========================================================="
    echo " "
    /BiO/Share/Tool/plink2 \
        --bfile {bfile} \
        --adjust \
        --ci 0.95 \
        --covar {covar} \
        --covar-name  {covar_name}\
        --linear hide-covar sex \
        --out {outprefix} \
        --pheno {pheno} \
        --pheno-name {pheno_name} \
        --covar-variance-standardize
    echo " "
    echo "=========================================================="
    echo "finished       : $(date)"
    echo "=========================================================="'''

    
    return cmd

# %%
bfile = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k/9000Korean_samples.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.excesshet_60.abhet_0.4.abhom_0.1.adsupport_0.9.kinship_3rd"
covar = "/BiO/Access/kyungwhan1998/genome/gwas/Resources/Data/KGP_age_sex_bmi_pc1_pc10_for_GWAS.txt"
covar_name = "age,age_sq,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
# covar_name = "age,age_sq,lymph,monopa,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
pheno = "/BiO/Access/kyungwhan1998/genome/gwas/Resources/Data/KGP_sbp_glucose_age_rdw_creat_bun_phenoage_calibrated_for_GWAS.txt"
outprefix = "/BiO/Access/kyungwhan1998/genome/gwas/Results/KGP_sbp_glucose_age_rdw_creat_bun_phenoage_calibrated_GWAS"
pheno_name = "phenoage_calibrated_advance_INT"

cmd = get_cmd(bfile, covar, covar_name, outprefix, pheno, pheno_name)
subprocess.run(cmd, shell=True)