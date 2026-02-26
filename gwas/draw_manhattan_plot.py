# %%
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# %%
path_assoc = "/BiO/Access/kyungwhan1998/genome/gwas/Results/KGP_creat_age_glucose_sbp_phenoage_calibrated_GWAS.phenoage_calibrated_advance_INT.glm.linear"
path_vep = "/BiO/Research/Korea10KGenome/Results/VariantAnnotation.by.VEP115/Korea10K_VariantAnnotation.FilterIndividualABHet2STD.VariantABHet0.4ABHom0.9ADSupp0.9.Biallelic.Geno_0.01.Mind_0.1.HWE_1e6.Het_3STD.Kinship_3rd.nonKorean.20251114.maf"
path_vep_gwas = "/BiO/Access/kyungwhan1998/genome/gwas/Results/KGP_creat_age_glucose_sbp_phenoage_calibrated_GWAS.phenoage_calibrated_advance_INT.glm.linear.annot"
gwas = pd.read_csv(path_assoc, sep="\t")
# %%
gwas_p_sorted = gwas[gwas["P"] < 1e-05].sort_values(by=["P"])

# %%
gwas_p_sorted["VAR"] = gwas_p_sorted["ID"] + "_" + gwas_p_sorted["REF"] + "_" + gwas_p_sorted["ALT"]

# %%
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

gwas_p_sorted["VAR_New"] = gwas_p_sorted["VAR"].apply(normalize_variant)

# %%
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

# %% 
new_header = ['Chromosome',
 'vcf_pos',
 'Reference_Allele',
 'Alternate_Allele',
 'dbSNP_RS',
 'Hugo_Symbol']

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

# %%
def manhattan_plot(
    df,
    chr_col='#CHROM',
    pos_col='POS',
    p_col='P',
    ax=None,
    sigline=5e-8,
    vep_df=None
):
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    
    list_annot_text = []
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(12,4))

    df = df.copy()
    df = df[df[chr_col].notnull()]

    # --- normalize chr column ---
    df[chr_col] = df[chr_col].astype(str).str.replace("chr", "", regex=False)

    chr_map = {'X': 23, 'Y': 24}
    df['CHRint'] = df[chr_col].replace(chr_map).astype(int)
    df['CHRstr'] = df['CHRint'].astype(str)

    # Lead SNPs per chromosome
    lead_idx = df.groupby('CHRint')[p_col].idxmin()
    df['is_lead'] = False
    df.loc[lead_idx, 'is_lead'] = True

    # Cumulative position
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

    # Colors
    unique_chrs = sorted(plot_df['CHRstr'].unique(), key=int)
    colors = ["0.3" if i%2==0 else "0.6" for i in range(len(unique_chrs))]
    color_map = {c: colors[i] for i, c in enumerate(unique_chrs)}
    plot_df['color'] = plot_df['CHRstr'].map(color_map)

    # -log10(P)
    plot_df['-log10p'] = -np.log10(plot_df[p_col].astype(float).clip(lower=1e-300))
    ax.scatter(plot_df['cumuBP'], plot_df['-log10p'], c=plot_df['color'], s=6)

    # -------------------------
    # Only annotate & plot lead SNPs
    # -------------------------
    lead_df = plot_df[plot_df['is_lead']]
    lead_df_suggestive = lead_df[lead_df[p_col] < sigline]

    # Red points ONLY for suggestive lead SNPs
    ax.scatter(
        lead_df_suggestive['cumuBP'],
        lead_df_suggestive['-log10p'],
        c='red',
        s=20
    )

    # Axes
    ax.set_xlabel("Chromosome")
    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels, rotation=90)
    ax.set_ylabel("-log10(P)")

    if sigline:
        ax.axhline(-np.log10(sigline), color='red', linestyle='--', linewidth=1)

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

        ann = merged[merged['Annotation'].notna()]

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
                fontsize=10,
                ha='center',
                color='red'
            )
            list_annot_text.append(annot_text)

    return ax, list_annot_text
# %%
import seaborn as sns
from adjustText import adjust_text

vep_gwas = pd.read_csv(path_vep_gwas, sep="\t")
vep_gwas["Annotation"] = vep_gwas["Hugo_Symbol"] + "\n" + "(" + vep_gwas["dbSNP_RS"] + ")"

ax, list_annot_text = manhattan_plot(gwas, sigline=1e-06, vep_df=vep_gwas)
ax.set_xlabel("Chromosome", fontsize=20)
ax.set_ylabel("-log10(P)", fontsize=18)
ax.set_xticklabels(ax.get_xticklabels(), fontsize=15)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=15)
# adjust_text(list_annot_text, arrowprops=dict(arrowstyle='->', color='red'))
sns.despine(ax=ax, top=True, right=True)
plt.show()
plt.close()

# %%
vep_gwas["VAR_New"] = vep_gwas["Chromosome"].astype(str) + "_" + vep_gwas["vcf_pos"].astype(str) + "_" + vep_gwas["Reference_Allele"] + "_" + vep_gwas["Alternate_Allele"]
 
gwas_vep_merged = gwas_p_sorted.merge(vep_gwas, how="inner", on="VAR_New")