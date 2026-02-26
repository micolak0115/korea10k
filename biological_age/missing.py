# %%
import numpy as np
import pandas as pd


# %%
class Impute():
    def __init__(self, imputing_data, imputed_data, list_target_genes):
        self.imputing_data = imputing_data
        self.imputed_data = imputed_data
        self.list_original_genes = list(self.imputed_data.columns)
        self.list_target_genes = list_target_genes

        if self.is_necessary_impute():
            self.list_imputed_genes = self.get_genes_need_imputation()
            print(f"Missing {len(self.list_imputed_genes)} genes: {', '.join(self.list_imputed_genes)}")

        else:
            print("Imputation Not Necessary!")

    def is_necessary_impute(self):
        num_need_impute = self.imputed_data.isnull().any().sum()
        if int(num_need_impute) > 0:

            return True

        else:
            
            return False

    def get_genes_need_imputation(self):
        list_col_imputing = list(self.imputing_data.columns)
        imputed_data = self.imputed_data.loc[:, self.imputed_data.notnull().any()]
        list_col_imputed = list(imputed_data.columns)
        
        list_imputed_genes = list(set(list_col_imputing).difference(set(list_col_imputed)))
        
        return list_imputed_genes

    def get_num_impute_rows(self):
        num_impute_rows = np.shape(self.imputing_data)[0]

        return num_impute_rows
    
    def get_num_impute_cols(self):
        num_impute_cols = len(self.list_imputed_genes)

        return num_impute_cols

    def zero_method(self):
        pass

    def mean_method(self):
        pass
    
    def knn_method(self):
        print("Imputing with the KNN method...")
        from sklearn.impute import KNNImputer
        merged_data = pd.concat([self.imputing_data, self.imputed_data], axis=0)
        imputer = KNNImputer(n_neighbors=min(10, np.shape(merged_data)[1]))
        impute_matrix = imputer.fit_transform(merged_data)
        df_impute = pd.DataFrame(impute_matrix, index=merged_data.index, columns=merged_data.columns)
        
        return df_impute

    def return_dataframe_filtered_imputed_samples_only(self, df_impute):
        list_sample_imputed_data = list(self.imputed_data.index)
        df_imputed = df_impute.loc[list_sample_imputed_data, :]
        
        return df_imputed