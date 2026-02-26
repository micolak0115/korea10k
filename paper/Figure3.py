#%%
import math
import re
from collections import Counter
from copy import deepcopy
from glob import glob

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import fdrcorrection

#%%
list_cpgrel_type = ["CG context Creation", "CG context Elimination", "Both Creation/Elimination", "Unrelated"]
dict_cpgrel_to_color = {
    "CG context Creation":"mediumblue",
    "CG context Elimination":"firebrick",
    "Both Creation/Elimination": "darkviolet",
    "Unrelated" : "darkgray"
}
dict_chisq_color = {
    "nonsig":"lightgray",
    "higher":"firebrick",
    "lower":"mediumblue"
}

#%%
population_compare = "Superpopulation_EUR"
window_size = 10000

path_var_gcount = "/BiO/Research/Korea10KGenome/Results/Plink.JointCall.to.hg38.with.AdapterTrimmedRead.for.BWA.mem.by.GATK.HaplotypeCaller.BoundaryMerged.RemoveABHetOutlier2STD.VQSR.PASS/Plink_Final/chromosome_merged.qc_filtered.kinship_filtered.nonKorean_filtered/Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.excesshet_60.abhet_0.4.abhom_0.1.adsupport_0.9.kinship_3rd.nonKorean.mac1.QC.gcount"
path_var_afreq = "/BiO/Research/Korea10KGenome/Results/Plink.JointCall.to.hg38.with.AdapterTrimmedRead.for.BWA.mem.by.GATK.HaplotypeCaller.BoundaryMerged.RemoveABHetOutlier2STD.VQSR.PASS/Plink_Final/chromosome_merged.qc_filtered.kinship_filtered.nonKorean_filtered/Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.excesshet_60.abhet_0.4.abhom_0.1.adsupport_0.9.kinship_3rd.nonKorean.mac1.QC.acount"
path_varnames_appeared = "/BiO/Research/Korea10KGenome/Workspace/Yoonsung/Results/Methyl_SNP_Comparison/SNP_related_to_CpG/QC_Filtered/nonKorean_mac1_Filtered/merged.CpG_appeared.Numeric_VCF.QC_Filtered.Filter_nonKorean.varnames"
path_varnames_disappeared = "/BiO/Research/Korea10KGenome/Workspace/Yoonsung/Results/Methyl_SNP_Comparison/SNP_related_to_CpG/QC_Filtered/nonKorean_mac1_Filtered/merged.CpG_disappeared.Numeric_VCF.QC_Filtered.Filter_nonKorean.varnames"

path_hg38 = "/BiO/Research/Korea10KGenome/Resources/Reference/chromosome/hg38.fa"
path_cpgs_appeared_format = "/BiO/Research/Korea10KGenome/Workspace/Yoonsung/Results/Methyl_SNP_Comparison/SNP_related_to_CpG/QC_Filtered/nonKorean_mac1_Filtered/Popular_Filtered/*.CpG_appeared.Numeric_VCF.QC_Filtered.Filter_nonKorean.Popular_Filtered_80.tsv"
path_cpgs_disappeared_format = "/BiO/Research/Korea10KGenome/Workspace/Yoonsung/Results/Methyl_SNP_Comparison/SNP_related_to_CpG/QC_Filtered/nonKorean_mac1_Filtered/Popular_Filtered/*.CpG_disappeared.Numeric_VCF.QC_Filtered.Filter_nonKorean.Popular_Filtered_80.tsv"

path_1kgp_comparison_result = "/BiO/Research/Korea10KGenome/Workspace/Yoonsung/Results/Methyl_SNP_Comparison/Compare_AF_with_Other_Population/Merged_chr.CpG_related.Summary.QC_Filtered.Filt_nonKorean.Popular_Filtered_80.window_10000.diff_10perc.AF_Compared_1KGP_Superpopulation_EUR.tsv"
#%%
def add_vartype_annotation(table):
    table_ = table.copy()
    table_snv = table_[np.logical_and(table_["REF"].apply(len) == 1, table_["ALT"].apply(len) == 1)]
    table_ins = table_[table_["ALT"].apply(len) > 1]
    table_del = table_[table_["REF"].apply(len) > 1]
    
    table_snv.loc[:, "VarType"] = "SNV"
    table_ins.loc[:, "VarType"] = "INS"
    table_del.loc[:, "VarType"] = "DEL"
    
    table_snv.loc[:, "Substitution"] = list(map(lambda val1, val2: f"{val1}>{val2}", table_snv["REF"], table_snv["ALT"]))
    
    table_ins.loc[:, "Length"] = table_ins["ALT"].apply(len) - 1
    table_del.loc[:, "Length"] = table_del["REF"].apply(len) - 1
    
    table_annot = pd.concat([table_snv, table_ins, table_del]).sort_values(by = ["#CHROM", "ID"]).reset_index(drop = True)
    return table_annot

def get_cpg_rel_status(is_excl, is_incl):
    if is_excl:
        if is_incl:
            return "Both Creation/Elimination"
        else:
            return "CG context Elimination" 
    elif is_incl:
        return "CG context Creation"
    else:
        return "Unrelated"

def read_fasta_as_dict(path_fasta):
    dict_fasta_to_lines = dict()

    header = None
    with open(path_fasta, 'r') as fr:
        for line in fr:
            if line.startswith('>'):
                header = line.strip('\n')[1:].split()[0]
                assert dict_fasta_to_lines.get(header) == None
                dict_fasta_to_lines[header] = list()
            else:
                dict_fasta_to_lines[header].append(line.strip('\n'))
    dict_fasta = dict()
    for name, lines in dict_fasta_to_lines.items():
        dict_fasta[name] = ''.join(lines)
    return dict_fasta

def read_cpg_related_snp_table_infoonly(path_table_format):
    list_paths = glob(path_table_format)
    
    list_tables = list(map(
        lambda path: pd.read_csv(path, sep = '\t', usecols = ["#CHROM", "POS", "ID", "REF", "ALT", "Bef_Nuc", "Next_Nuc", "Pos_C", "Pos_G"]),
        list_paths
    ))
    table_all = pd.concat(list_tables)
    return table_all

def add_ncg_col(table, col_check, col_add):
    re_cg = re.compile("CG")
    table[col_add] = table[col_check].apply(lambda val: len(re_cg.findall(val)))
    return table

def get_cg_index_from_fasta(dict_fasta, list_chrnames_check):
    re_cg = re.compile("CG")
    
    dict_chrname_to_cg_index = dict()
    for chrname in list_chrnames_check:
        list_match_result = list(re_cg.finditer(dict_fasta[chrname]))
        list_ind_cgs = list(map(lambda result: result.span()[0], list_match_result))
        
        dict_chrname_to_cg_index[chrname] = list_ind_cgs
    return dict_chrname_to_cg_index

def get_count_per_window_from_index(list_index, window_size):
    dict_window_ind_to_cnt = dict()
    for val in list_index:
        window_index = math.floor(val / window_size)
        
        if dict_window_ind_to_cnt.get(window_index) == None:
            dict_window_ind_to_cnt[window_index] = 0
        dict_window_ind_to_cnt[window_index] += 1
    return dict_window_ind_to_cnt

def get_variation_per_window_from_table(table, window_size, col_chr, col_pos, col_variation):
    dict_chr_to_ncg_diff_per_window = dict()
    for chrname in table[col_chr].unique():
        dict_chr_to_ncg_diff_per_window[chrname] = dict()
        
        table_chr = table[table[col_chr] == chrname]
        
        list_window_index = table_chr[col_pos].apply(lambda val: math.floor(val / window_size))
        
        for window_index, n_diff in zip(list_window_index, table_chr[col_variation]):
            if dict_chr_to_ncg_diff_per_window[chrname].get(window_index) == None:
                dict_chr_to_ncg_diff_per_window[chrname][window_index] = 0
            dict_chr_to_ncg_diff_per_window[chrname][window_index] += n_diff
    return dict_chr_to_ncg_diff_per_window

def run_chisquare_test(altct1, obsct1, altct2, obsct2):
    refct1 = obsct1 - altct1
    refct2 = obsct2 - altct2
    
    if refct1 == 0 and refct2 == 0:
        p = 1
    else:
        try:
            chi, p, dof, expected = chi2_contingency(np.array([[altct1, refct1], [altct2, refct2]]))
        except Exception as e:
            p = None
    return p

def get_color(fdr, diff):
    if fdr > 0.05:
        return dict_chisq_color["nonsig"]
    elif diff > 0:
        return dict_chisq_color["higher"]
    else:
        return dict_chisq_color["lower"]
#%%
table_gcount = pd.read_csv(path_var_gcount, delim_whitespace = True)
table_afreq = pd.read_csv(path_var_afreq, delim_whitespace = True)
table_afreq["ALT_FREQS"] = table_afreq["ALT_CTS"] / table_afreq["OBS_CT"]
list_varnames_appeared = pd.read_csv(path_varnames_appeared, names = ["varname"])["varname"].to_list()
list_varnames_disappeared = pd.read_csv(path_varnames_disappeared, names = ["varname"])["varname"].to_list()

dict_hg38 = read_fasta_as_dict(path_hg38)
table_appeared = read_cpg_related_snp_table_infoonly(path_cpgs_appeared_format)
table_disappeared = read_cpg_related_snp_table_infoonly(path_cpgs_disappeared_format)

table_snps_check_fdrcalc = pd.read_csv(path_1kgp_comparison_result, sep = '\t')
#%%
table_afreq_vartype_annot = add_vartype_annotation(table_afreq)

set_varnames_appeared = set(list_varnames_appeared)
set_varnames_disappeared = set(list_varnames_disappeared)

table_afreq_vartype_annot["CpG_Elimination"] = table_afreq_vartype_annot["ID"].apply(lambda val: val in set_varnames_disappeared)
table_afreq_vartype_annot["CpG_Creation"] = table_afreq_vartype_annot["ID"].apply(lambda val: val in set_varnames_appeared)

table_afreq_vartype_annot["CpG_Related_Status"] = list(map(lambda excl, incl: get_cpg_rel_status(excl, incl), table_afreq_vartype_annot["CpG_Elimination"], table_afreq_vartype_annot["CpG_Creation"]))


table_gcount["NONREF_CT"] = table_gcount["HET_REF_ALT_CTS"] + table_gcount["TWO_ALT_GENO_CTS"]
table_gcount["NONREF_RATIO"] = table_gcount["NONREF_CT"] / (table_gcount["NONREF_CT"] + table_gcount["HOM_REF_CT"])

set_all_varnames_related_to_cg = set_varnames_appeared|set_varnames_disappeared
table_gcount["CpG_Relation"] = table_gcount["ID"].apply(lambda val: val in set_all_varnames_related_to_cg)
table_gcount_cgv = table_gcount[table_gcount["CpG_Relation"]]

#%%
table_appeared["REF_Context"] = table_appeared["Bef_Nuc"] + table_appeared["REF"] + table_appeared["Next_Nuc"]
table_appeared["ALT_Context"] = table_appeared["Bef_Nuc"] + table_appeared["ALT"] + table_appeared["Next_Nuc"]

table_disappeared["REF_Context"] = table_disappeared["Bef_Nuc"] + table_disappeared["REF"] + table_disappeared["Next_Nuc"]
table_disappeared["ALT_Context"] = table_disappeared["Bef_Nuc"] + table_disappeared["ALT"] + table_disappeared["Next_Nuc"]

table_appeared = add_ncg_col(table_appeared, "REF_Context", "REF_N_CG_Context")
table_appeared = add_ncg_col(table_appeared, "ALT_Context", "ALT_N_CG_Context")

table_disappeared = add_ncg_col(table_disappeared, "REF_Context", "REF_N_CG_Context")
table_disappeared = add_ncg_col(table_disappeared, "ALT_Context", "ALT_N_CG_Context")

table_appeared["N_Created"] = table_appeared["ALT_N_CG_Context"] - table_appeared["REF_N_CG_Context"]
table_disappeared["N_Created"] = table_disappeared["ALT_N_CG_Context"] - table_disappeared["REF_N_CG_Context"]

table_appeared["ID"] = table_appeared["#CHROM"] + ':' + table_appeared["POS"].astype(str)
table_disappeared["ID"] = table_disappeared["#CHROM"] + ':' + table_disappeared["POS"].astype(str)

table_cpg_rel_snp_total = pd.concat([table_appeared, table_disappeared]).drop_duplicates(subset = ["ID"]).sort_values(by = ["#CHROM", "POS"]).reset_index(drop = True)

list_chrnames_check = list(map(lambda val: f"chr{val}", range(1, 23)))
dict_hg38_chrname_to_cg_index = get_cg_index_from_fasta(dict_hg38, list_chrnames_check)

dict_hg38_len = dict()
for chrname in list_chrnames_check:
    dict_hg38_len[chrname] = len(dict_hg38[chrname])

dict_hg38_chrname_to_ncg_per_window = dict()
for chrname in list_chrnames_check:
    dict_hg38_chrname_to_ncg_per_window[chrname] = get_count_per_window_from_index(dict_hg38_chrname_to_cg_index[chrname], window_size)

dict_chrname_to_ncg_diff_per_window = get_variation_per_window_from_table(table_cpg_rel_snp_total, window_size, "#CHROM", "POS", "N_Created")

dict_chrname_to_ratio_diff_per_window = dict()
for chrname in list_chrnames_check:
    dict_chrname_to_ratio_diff_per_window[chrname] = dict()
    
    max_window = math.floor(dict_hg38_len[chrname] / window_size)
    
    for ind_window in range(max_window+1):
        n_cg_hg38 = dict_hg38_chrname_to_ncg_per_window[chrname].get(ind_window, 0)
        n_cg_diff = dict_chrname_to_ncg_diff_per_window[chrname].get(ind_window, 0)
        
        if n_cg_hg38 == 0:
            if n_cg_diff == 0:
                dict_chrname_to_ratio_diff_per_window[chrname][ind_window] = np.nan
            else:
                dict_chrname_to_ratio_diff_per_window[chrname][ind_window] = np.inf
                print(f"{chrname} {ind_window} : inf")
        else:
            dict_chrname_to_ratio_diff_per_window[chrname][ind_window] = n_cg_diff / n_cg_hg38

#%%
table_snps_check_fdrcalc["AF_Diff"] = table_snps_check_fdrcalc["AF_Korea"] - table_snps_check_fdrcalc["AF_1KGP"]

table_snps_check_fdrcalc_region_increased = table_snps_check_fdrcalc[table_snps_check_fdrcalc["Ratio_Diff"] > 0]
table_snps_check_fdrcalc_region_decreased = table_snps_check_fdrcalc[table_snps_check_fdrcalc["Ratio_Diff"] < 0]

table_snps_check_fdrcalc_region_increased_positive = table_snps_check_fdrcalc_region_increased[table_snps_check_fdrcalc_region_increased["N_Created"] > 0]
table_snps_check_fdrcalc_region_decreased_negative = table_snps_check_fdrcalc_region_decreased[table_snps_check_fdrcalc_region_decreased["N_Created"] < 0]

table_draw_ax1 = table_snps_check_fdrcalc_region_increased_positive.copy()
table_draw_ax2 = table_snps_check_fdrcalc_region_decreased_negative.copy()

table_draw_ax1["Log10FDR"] = table_draw_ax1["ChiSquared_FDR"].apply(lambda val: -np.log10(val + 1e-320))
table_draw_ax2["Log10FDR"] = table_draw_ax2["ChiSquared_FDR"].apply(lambda val: -np.log10(val + 1e-320))

table_draw_ax1["color"] = table_draw_ax1.apply(lambda row: get_color(row["ChiSquared_FDR"], row["AF_Diff"]), axis = 1)
table_draw_ax2["color"] = table_draw_ax2.apply(lambda row: get_color(row["ChiSquared_FDR"], row["AF_Diff"]), axis = 1)
#%%
list_vartype = ["SNV", "INS", "DEL"]

dict_vartype_to_cpgrel_count = dict()
for vartype in list_vartype:
    dict_vartype_to_cpgrel_count[vartype] = Counter(table_afreq_vartype_annot[table_afreq_vartype_annot["VarType"] == vartype]["CpG_Related_Status"])

list_total_ratio = list()
list_total_diff = list()
for chrname in dict_chrname_to_ratio_diff_per_window.keys():
    list_index = list(dict_chrname_to_ratio_diff_per_window[chrname].keys())
    list_ratios = list(map(dict_chrname_to_ratio_diff_per_window[chrname].__getitem__, list_index))
    list_cnts = list(map(lambda ind: dict_chrname_to_ncg_diff_per_window[chrname].get(ind, np.nan), list_index))
    
    list_zip_ratio_cnt = list(zip(list_ratios, list_cnts))
    list_zip_ratio_cnt_nonna = list(filter(lambda val: pd.notna(val[0]), list_zip_ratio_cnt))
    
    list_total_ratio.extend(list(map(lambda val: val[0], list_zip_ratio_cnt_nonna)))
    list_total_diff.extend(list(map(lambda val: val[1], list_zip_ratio_cnt_nonna)))

#%%
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams["font.size"] = 17
fontsize_label = 25

fig = plt.figure(figsize=(10, 15))
row = 150
col = 101
gsfig = GridSpec(
    row, col, 
    left=0, right=1, bottom=0,
    top=1, wspace=1, hspace=1)

layer1_width = 27
layer1_wspace = 10

gs_cgv_ac = gsfig[0:40, 0:27]
ax_cgv_ac = fig.add_subplot(gs_cgv_ac)

gs_cgv_ac_ratio = gsfig[0:40, 35:35+27]
ax_cgv_ac_ratio = fig.add_subplot(gs_cgv_ac_ratio)

gs_popular_cgv = gsfig[50:85, 0:45]
ax_popular_cgv = fig.add_subplot(gs_popular_cgv)

gs_diff2ratio = gsfig[50:85, col-45:col]
ax_diff2ratio = fig.add_subplot(gs_diff2ratio)

gs_cgv_c = gsfig[95:150, 0:45] 
ax_cgv_c = fig.add_subplot(gs_cgv_c)

gs_cgv_e = gsfig[95:150, col-45:col]
ax_cgv_e = fig.add_subplot(gs_cgv_e)

label_x = -0.1
label_y = 0.02


############# CGV on Korea10K
bottom = np.array([0] * len(list_vartype))
for cpgreltype in list_cpgrel_type:
    list_cnt = np.array(list(map(lambda vartype: dict_vartype_to_cpgrel_count[vartype][cpgreltype], list_vartype)))
    
    ax_cgv_ac.bar(
        range(len(list_vartype)),
        list_cnt,
        color = dict_cpgrel_to_color[cpgreltype],
        bottom = bottom,
        # label = cpgreltype,
        width = 0.8
    )
    print(list_cnt)
    bottom += list_cnt
ax_cgv_ac.set_xticks(range(len(list_vartype)), list_vartype)
ax_cgv_ac.set_ylabel("# Variants")
ax_cgv_ac.set_xlim(-0.7, 2.7)
ax_cgv_ac.text(
    label_x, 1+label_y, 'A', fontsize = fontsize_label, fontweight = "bold", transform = ax_cgv_ac.transAxes, ha = "right", va = "bottom"
)

bottom = np.array([0.0] * len(list_vartype))
for cpgreltype in list_cpgrel_type:
    list_cnt = np.array(list(map(lambda vartype: dict_vartype_to_cpgrel_count[vartype][cpgreltype], list_vartype)))
    list_ratio = list_cnt / np.array(list(map(lambda vartype: sum(dict_vartype_to_cpgrel_count[vartype].values()), list_vartype)))
    
    ax_cgv_ac_ratio.bar(
        range(len(list_vartype)),
        list_ratio,
        color = dict_cpgrel_to_color[cpgreltype],
        bottom = bottom,
        label = cpgreltype,
        width = 0.8
    )
    bottom += list_ratio
ax_cgv_ac_ratio.set_xticks(range(len(list_vartype)), list_vartype)
ax_cgv_ac_ratio.set_ylabel("Ratio of Variants")
ax_cgv_ac_ratio.set_xlim(-0.7, 2.7)
ax_cgv_ac_ratio.text(
    label_x * ((gs_cgv_ac.colspan.stop-gs_cgv_ac.colspan.start) / (gs_cgv_ac_ratio.colspan.stop-gs_cgv_ac_ratio.colspan.start)), 1+label_y*((gs_cgv_ac.rowspan.stop-gs_cgv_ac.rowspan.start) / (gs_cgv_ac_ratio.rowspan.stop-gs_cgv_ac_ratio.rowspan.start)), 'B', fontsize = fontsize_label, fontweight = "bold", transform = ax_cgv_ac_ratio.transAxes, ha = "right", va = "bottom"
)

ax_cgv_ac_ratio.legend(loc = "center left", bbox_to_anchor = (1.01, 0.5), ncol = 1)


################ Popular variant from CGV
ax_popular_cgv.hist(table_gcount[table_gcount["CpG_Relation"]]["NONREF_RATIO"], bins = np.linspace(0, 1, 101), color = "gray")
ax_popular_cgv.set_yscale("log")
ax_popular_cgv.set_xlabel("Ratio of Non-REF Genotype")
ax_popular_cgv.set_ylabel("# Variant")
ax_popular_cgv.axvline(0.8, linewidth = 2, color = "firebrick", linestyle = "--")
ax_popular_cgv.text(
    label_x * ((gs_cgv_ac.colspan.stop-gs_cgv_ac.colspan.start) / (gs_popular_cgv.colspan.stop-gs_popular_cgv.colspan.start)), 1+label_y*((gs_cgv_ac.rowspan.stop-gs_cgv_ac.rowspan.start) / (gs_popular_cgv.rowspan.stop-gs_popular_cgv.rowspan.start)), 'C', fontsize = fontsize_label, fontweight = "bold", transform = ax_popular_cgv.transAxes, ha = "right", va = "bottom"
)

################ CG diff to Ratio
ax_diff2ratio.hist(list_total_ratio, bins = np.linspace(-0.3, 0.3, 49), color = "gray")
ax_diff2ratio.set_xlabel("Ratio of # CG Difference")
ax_diff2ratio.set_ylabel("# Windows")
# ax_diff2ratio.set_yscale("log")
# ax_diff2ratio.scatter(list_total_diff, list_total_ratio, s = 3, color = "gray")
# ax_diff2ratio.set_xlabel("# CG Difference")
# ax_diff2ratio.set_ylabel("Ratio of # CG Difference")
ax_diff2ratio.text(
    label_x * ((gs_cgv_ac.colspan.stop-gs_cgv_ac.colspan.start) / (gs_diff2ratio.colspan.stop-gs_diff2ratio.colspan.start)), 1+label_y*((gs_cgv_ac.rowspan.stop-gs_cgv_ac.rowspan.start) / (gs_diff2ratio.rowspan.stop-gs_diff2ratio.rowspan.start)), 'D', fontsize = fontsize_label, fontweight = "bold", transform = ax_diff2ratio.transAxes, ha = "right", va = "bottom"
)

################ 1KGP Diff
ax_cgv_c.scatter(
    table_draw_ax1["AF_Diff"],
    table_draw_ax1["Log10FDR"],
    s = 3,
    color = table_draw_ax1["color"],
    zorder = 1
)
ax_cgv_c.axhline(-np.log10(0.05), linewidth = 1, linestyle = "--", color = "gray", zorder = 3)
ax_cgv_c.set_xlabel("AF Difference")
ax_cgv_c.set_ylabel("-Log10(FDR)")
ax_cgv_c.set_xlim(-1, 1)
ax_cgv_c.text(
    label_x * ((gs_cgv_ac.colspan.stop-gs_cgv_ac.colspan.start) / (gs_cgv_c.colspan.stop-gs_cgv_c.colspan.start)), 1+label_y*((gs_cgv_ac.rowspan.stop-gs_cgv_ac.rowspan.start) / (gs_cgv_c.rowspan.stop-gs_cgv_c.rowspan.start)), 'E', fontsize = fontsize_label, fontweight = "bold", transform = ax_cgv_c.transAxes, ha = "right", va = "bottom"
)

ax_cgv_e.scatter(
    table_draw_ax2["AF_Diff"],
    table_draw_ax2["Log10FDR"],
    s = 3,
    color = table_draw_ax2["color"],
    zorder = 1
)
ax_cgv_e.axhline(-np.log10(0.05), linewidth = 1, linestyle = "--", color = "gray", zorder = 3)
ax_cgv_e.set_xlabel("AF Difference")
ax_cgv_e.set_ylabel("-Log10(FDR)")
ax_cgv_e.set_xlim(-1, 1)
ax_cgv_e.text(
    label_x * ((gs_cgv_ac.colspan.stop-gs_cgv_ac.colspan.start) / (gs_cgv_e.colspan.stop-gs_cgv_e.colspan.start)), 1+label_y*((gs_cgv_ac.rowspan.stop-gs_cgv_ac.rowspan.start) / (gs_cgv_e.rowspan.stop-gs_cgv_e.rowspan.start)), 'F', fontsize = fontsize_label, fontweight = "bold", transform = ax_cgv_e.transAxes, ha = "right", va = "bottom"
)

plt.show()