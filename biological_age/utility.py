# %%
import gzip
import os
import pickle

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler


# %%
class Utility():
    def __init__(self):
        self.gene_tag = "ENSG"
        self.remove_genes = ["PAR_Y", "CTC-338M12.4"]
        self.verbose = False
        self.returnX = False
        self.returnY = True
        self.sep = "\t"

    def read_metadata(self, inmeta, colsample):
        df_meta = pd.read_csv(inmeta, sep="\t")
        first_colheader = list(df_meta.columns)[0]
        if first_colheader != colsample:
            df_meta = df_meta.rename(columns={first_colheader: colsample})
        
        return df_meta

    def save_metadata(self, df_metadata, path_output):
        df_metadata.to_csv(path_output, sep="\t", index=False)

        return None

    def select_samples_meta(self, df_meta, list_select_samples, colsample):
        df_meta_select = df_meta[df_meta[colsample].isin(list_select_samples)]

        return df_meta_select  

    def drop_samples_meta(self, df_meta, list_drop_samples, colsample):
        df_meta_drop = df_meta[~df_meta[colsample].isin(list_drop_samples)]

        return df_meta_drop

    def read_expmtx(self, expmtx):
        df_exp = pd.read_csv(expmtx, sep="\t")
        
        return df_exp

    def drop_samples_exp(self, df_exp, list_drop_samples):
        df_exp_drop = df_exp.drop(columns=list_drop_samples)

        return df_exp_drop

    def select_samples_df(self, df_exp, list_select_samples, colgene):
        df_exp_indexed = df_exp.set_index(colgene)
        df_exp_select = df_exp_indexed[list_select_samples]
        df_exp_select_reidx = df_exp_select.reset_index(drop=False)

        return df_exp_select_reidx

    def save_expmtx(self, df_save, path_output):
        df_save.to_csv(path_output, sep="\t", index=False)

        return None

    def get_list_drop_samples(self, file_drop_list):
        with open(file_drop_list, mode="r") as fr:
            list_drop_samples = fr.readlines()
        list_drop_samples = list(map(lambda x: x.rstrip("\n"), list_drop_samples))
        
        return list_drop_samples

    def get_file_suffix_age_group(self, list_age_group_selected):
        list_age_group_selected = list(map(str, list_age_group_selected))
        file_suffix = '_'.join(list_age_group_selected)
        
        return file_suffix
        
    def add_df_sample_age_group(self, df_meta):
        df_meta["Sample_Age_Group"] = df_meta["Sample_Age"].map(lambda x: int(str(x)[0] + "0"))
        
        return df_meta

    def remove_df_certain_age_group(self, df_meta, list_excl_age_group):
        df_meta_filtered_out = df_meta[df_meta["Sample_Age_Group"].isin(list_excl_age_group)]
        
        return df_meta_filtered_out

    def get_sampleid_removed_age_group(self, df_meta_filtered_out, colsample):
        list_samples_excl = df_meta_filtered_out[colsample].to_list()
        
        return list_samples_excl

    def filter_out_df_excl_samples(self, df_meta_ori, list_samples_excl, colsample):
        df_meta_filtered_in = df_meta_ori[~df_meta_ori[colsample].isin(list_samples_excl)]
        
        return df_meta_filtered_in

    def get_list_column_val(self, df_meta_filtered_in, column="Sample_Age_Group"):
        list_column_val = df_meta_filtered_in[column].to_list()
        
        return list_column_val

    def get_dictionary_sample_age_group(self, list_samples, list_age_group):
        dict_samples_age_group = dict(zip(list_samples, list_age_group))
        
        return dict_samples_age_group

    def get_randomly_selected_samples_certain_age_group(self, dict_samples_age_group, list_age_group_selected, num_samples_selected=100, random_seed=1):
        dict_samples_age_group_selected = {sample: age_group for sample, age_group in dict_samples_age_group.items() if int(age_group) in list_age_group_selected}
        num_samples_age_group_selected = len(dict_samples_age_group_selected.keys())
        if num_samples_selected > num_samples_age_group_selected:
            num_samples_selected = num_samples_age_group_selected
            print(f"Warning: Less samples selected - only {num_samples_age_group_selected} samples in the age group(s)!")
            print(f"Proceeding with non-random selection...")
            selected_samples = list(dict_samples_age_group_selected.keys())
        else:
            np.random.seed(random_seed)
            selected_samples = list(np.random.choice(list(dict_samples_age_group_selected.keys()), size=num_samples_selected, replace=False))
        
        return selected_samples

    def get_list_sample_id(self, inexp: str, colgene: str) -> str:
        list_sample_id = list()
        with open(inexp, mode="r") as fr:
            sampleid = fr.readline().rstrip("\n").split("\t")
            list_sample_id.extend(sampleid)
        if colgene in list_sample_id:
            list_sample_id.remove(colgene)

        return list_sample_id

    def get_dictionary_sample_age(self, inmeta: str, colsample: str, colage: str) -> str:
        dict_sample_age = dict()
        with open(inmeta, mode="r") as fr_meta:
            header = fr_meta.readline()
            idx_sample = header.index(colsample)
            idx_age = header.index(colage)
            for line in fr_meta:
                record = line.rstrip("\n").split("\t")
                sampleid = record[idx_sample]
                age = record[idx_age]

                if len(age.split("-")) > 0:
                    age = (int(age.split("-")[0]) + int(age.split("-")[1]) + 1) / 2

                dict_sample_age[sampleid] = age
        
        return dict_sample_age

    def get_dictionary_gene_expression(self, inexp: str) -> dict:
        dict_gene_exp = dict()
        with open(inexp, mode="r") as fr:
            header = fr.readline()
            for line in fr:
                record = line.rstrip("\n").split("\t")
                gene = record[0]
                expval = record[1:]
                dict_gene_exp[gene] = expval
        
        return dict_gene_exp

    def get_dictionary_age_expression(self, dict_gene_exp: dict, dict_sample_age: dict):
        dict_age_exp = dict()
        for sampleexp, exp in dict_gene_exp.items():
            for samplemeta, age in dict_sample_age.items():
                if sampleexp == samplemeta:
                    dict_age_exp[age] = list()
                    dict_age_exp[age].extend(exp)
        
        return dict_age_exp

    def get_list_median_expression_per_gene(self, dict_gene_exp_filt: dict):
        list_medianexp = list()
        for _, list_exp in dict_gene_exp_filt.items():
            list_exp_float = list(map(float, list_exp))
            medianexp = np.median(list_exp_float)
            list_medianexp.append(medianexp) 
        
        return list_medianexp

    def make_feature_table(self, dict_gene_exp_filt: dict, list_sample_id: list, colsample:str) -> pd.DataFrame:
        df_feature = pd.DataFrame.from_records(dict_gene_exp_filt)
        df_feature.index = list_sample_id
        df_feature_reset_idx = df_feature.reset_index(drop=False)
        df_feature_reset_idx_rename = df_feature_reset_idx.rename(columns={"index":colsample})
        
        return df_feature_reset_idx_rename

    def merge_feature_table_with_metadata(self, df_feature: pd.DataFrame, df_meta: pd.DataFrame, merge_how: str, merge_on: str) -> pd.DataFrame:
        df_feature_meta_merged = pd.merge(df_feature, df_meta, how=merge_how, on=merge_on)

        return df_feature_meta_merged

    def save_merged_feature_table(self, df_save: pd.DataFrame, path_feature:str) -> pd.DataFrame:
        df_save.to_csv(path_feature, sep="\t", index=False)

    def get_stratify_bins(self, y, nbins):
        step = (max(y) - min(y)) / nbins
        print(nbins)
        print(step)
        bins = np.arange(min(y), max(y)+step/2, step)
        bins[0] -= step/2
        bins[-1] += step/2
        print(max(y))
        print(min(y))
        print(bins)
        stratify = np.digitize(y, bins=bins, right=True)
        
        return stratify

    def stratify_split_train_test_data(self, X: pd.DataFrame, y: pd.DataFrame, stratify=None, train_ratio=0.70, random_state=1) -> pd.DataFrame:
        X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=stratify, test_size= 1 - train_ratio, random_state=random_state)
        
        util = Utility()
        X_train_reidx = util.reset_indices(X_train)
        X_test_reidx = util.reset_indices(X_test)
        y_train_reidx = util.reset_indices(y_train)
        y_test_reidx = util.reset_indices(y_test)
        
        return X_train_reidx, X_test_reidx, y_train_reidx, y_test_reidx

    def dump_pickle(self, data: str, dir_out: str, filename: str) ->  str:
        outfilename = os.path.join(dir_out, filename)
        os.makedirs(os.path.dirname(outfilename), exist_ok=True)
        with gzip.open(outfilename, 'wb') as f:
            pickle.dump(data, f)

        return None

    def load_pickle(self, dir_out: str, filename: str) -> str:
        loadfilename = os.path.join(dir_out, filename)
        with gzip.open(loadfilename,'rb') as f:
            data = pickle.load(f)
        
        return data
    
    def reset_indices(self, df: pd.DataFrame) -> pd.DataFrame:
        df_reidx = df.reset_index(drop=True)

        return df_reidx

    def get_end_idx_gene_list(self, df_feature: pd.DataFrame, gene_tag: str) -> int:
        list_header = list(df_feature.columns)
        list_gene_idx = list()
        for idx, column in enumerate(list_header):
            if str(column).startswith(gene_tag):
                list_gene_idx.append(idx)
        end_idx_gene = int(max(list_gene_idx))

        return end_idx_gene

    def lookup_available_variables(self, df_dataset: pd.DataFrame) -> list:
        list_dataset = df_dataset.columns.to_list()
        
        return list_dataset

    def overview_variables(self, X: pd.DataFrame, ys: pd.DataFrame) -> list:
        list_x = self.lookup_available_variables(X)
        print(f"Xs: {list_x}")

        list_y = self.lookup_available_variables(ys)
        print(f"Ys: {list_y}")

        if self.returnX:
            return list_x

        if self.returnY:
            return list_y
        
        if (self.returnX and self.returnY):
            return list_x, list_y

    def get_list_selected_genes(self, file_linear_reg):
        list_selected_genes = list()
        with open(file_linear_reg, mode="r") as fr:
            _skiprow = fr.readline()
            for line in fr:
                record = line.rstrip("\n").split("\t")
                gene = record[0]
                list_selected_genes.append(gene)

        return list_selected_genes

    def save_selected_genes(self, dir_out, list_selected_genes):
        fileingenename = os.path.join(dir_out, "features_ingenename.txt")
        with open(fileingenename, mode="w") as fw:
            for gene in list_selected_genes:
                fw.write(gene + "\n")
    
    def get_feature_table(self, feature_table) -> pd.DataFrame:
        df_feat = pd.read_csv(feature_table, sep=self.sep, low_memory=False)
        end_idx_gene = self.get_end_idx_gene_list(df_feat, self.gene_tag)
        X = df_feat.iloc[:, :(end_idx_gene+1)]
        for remove_tag in self.remove_genes:
            list_filt_colnames = list(filter(lambda x: remove_tag not in x, list(X.columns[1:])))
        X_filtered = X[[X.columns[0]] + list_filt_colnames]
        list_new_colnames = [X_filtered.columns[0]] + list(map(lambda x: x.split(".")[0] + "_" + "_".join(x.split("_")[1:]) if len(x.split("."))>1 else x, list_filt_colnames))
        X_filtered.columns = list_new_colnames
        ys = df_feat.iloc[:, (end_idx_gene+1):]

        if self.verbose:
            self.overview_variables(X, ys)

        return X_filtered, ys

    def filter_input_data(self, X, ys, col_y):
        list_idx_null = list(np.where(~pd.isnull(ys[col_y]))[0])
        y = ys[col_y].iloc[list_idx_null]
        X_ = X.iloc[list_idx_null, :]
        print(f"'{col_y}' selected")

        return X_, y

    def get_list_gene_name(self, X):
        list_gene_name = list(X.iloc[:, 1:].columns)

        return list_gene_name
    
    def get_target_genes(self, path_target):
        with open(path_target, mode="r") as fr:
            list_target_genes = fr.readlines()
        
        list_target_genes = list(map(lambda x: x.rstrip("\n"), list_target_genes))

        return list_target_genes
    
    def select_dataset_target_genes(self, X_set_idx, list_target_genes, colsample):
        for gene in list_target_genes:
            if gene not in X_set_idx.columns:
                X_set_idx[gene] = np.NaN
        X_target = X_set_idx[list_target_genes]
            
        return X_target
    
    def standardize_expression(self, X_train_set_idx, X_test_set_idx, colsample):        
        scaler = StandardScaler().fit(X_train_set_idx)

        list_train_samples = list(X_train_set_idx.index)
        list_train_genes = list(X_train_set_idx.columns)

        list_test_samples = list(X_test_set_idx.index)
        list_test_genes = list(X_test_set_idx.columns)

        # Apply standardization
        X_train_set_idx_std = pd.DataFrame(scaler.transform(X_train_set_idx), 
                                        index=list_train_samples, 
                                        columns=list_train_genes)

        X_test_set_idx_std = pd.DataFrame(scaler.transform(X_test_set_idx), 
                                        index=list_test_samples, 
                                        columns=list_test_genes)

        # Reset index and rename the index column
        X_train_reset_idx_std = X_train_set_idx_std.reset_index(drop=False).rename(columns={"index": colsample})
        X_test_reset_idx_std = X_test_set_idx_std.reset_index(drop=False).rename(columns={"index": colsample})

        # Return the standardized data
        return X_train_reset_idx_std, X_test_reset_idx_std
    

    def save_input_dataset(self, X_, y_, outfilename, outdir, **kwargs):
        testing_dataset = pd.concat([X_, y_], axis=1)
        testing_dataset.to_csv(os.path.join(outdir, outfilename), sep=self.sep, **kwargs)

        return testing_dataset
# %%
