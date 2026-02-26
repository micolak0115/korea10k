#%%
## Courtesy of Yoonsung (2025-01-22)
## /BiO/Access/yoonsung/SpyWare/PhenotypicAge/Resources/Scripts/run_lifelines_cv.py

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from lifelines import CoxPHFitter
from matplotlib import pyplot as plt
from sklearn.model_selection import KFold


# %%
def select_best_penalizer(table, col_duration, col_event, penalizer_range, l1_ratio, cv = 5, random_state = 0, n_jobs = 10):
    list_dict_param = get_parameter_combinations(penalizer_range, l1_ratio)
    
    with Parallel(n_jobs = n_jobs) as parallel:
        list_of_list_dev = parallel(
            delayed(get_deviances_with_specific_parameter_cv)(
                table, col_duration, col_event, params, cv, random_state
            ) for params in list_dict_param
        )
    
    ind_min, ind_min_std, ind_min_half_std = find_best_penalizer_from_deviances(list_of_list_dev)
    
    dict_result = {
        "Partial_LLH_deviance" : list_of_list_dev,
        "Penalizers" : penalizer_range,
        "l1_ratio" : l1_ratio,
        "Index_min_deviance" : ind_min,
        "Index_min_plus_std_deviance" : ind_min_std,
        "Index_min_plus_half_std_deviance": ind_min_half_std
    }
    
    return dict_result    
    
def get_parameter_combinations(penalizer_range, l1_ratio):
    list_dict_param = list()
    for pen in penalizer_range:
        list_dict_param.append({
            "penalizer":pen,
            "l1_ratio":l1_ratio
        })
    return list_dict_param
    
def get_deviances_with_specific_parameter_cv(table, col_duration, col_event, params, cv, random_state):
    kfold_split = KFold(n_splits=cv, shuffle = True, random_state=random_state)
    list_train_test_zip = list()
    for index_train, index_test in kfold_split.split(table):
        table_train_i = table.loc[index_train, :]
        table_test_i = table.loc[index_test, :]
        list_train_test_zip.append((table_train_i, table_test_i))
    
    list_deviance = list()
    for table_train, table_test in list_train_test_zip:
        partial_loglikelihood_deviance = get_testset_deviance_with_specific_paramter(table_train, table_test, col_duration, col_event, params)
        list_deviance.append(partial_loglikelihood_deviance)
    return list_deviance    
    
def get_testset_deviance_with_specific_paramter(table_train, table_test, col_duration, col_event, params):
    cphmodel = CoxPHFitter(**params)
    cphmodel.fit(table_train, duration_col=col_duration, event_col=col_event)
    
    partial_loglikelihood = cphmodel.score(table_test)
    partial_loglikelihood_deviance = -2 * partial_loglikelihood
    
    return partial_loglikelihood_deviance
    
def find_best_penalizer_from_deviances(list_of_list_dev):
    mean_dev_per_penalizer = np.mean(list_of_list_dev, axis = 1)
    ind_min_mean_dev = np.argmin(mean_dev_per_penalizer)
    minimum_mean_dev = min(mean_dev_per_penalizer)
    
    std_of_min_dev = np.std(list_of_list_dev[ind_min_mean_dev])
    dev_allow_range = minimum_mean_dev+std_of_min_dev
    list_ind_less_than_1std_from_min = list(filter(lambda ind: mean_dev_per_penalizer[ind]<dev_allow_range, range(len(list_of_list_dev))))
    ind_maximum_dev_from_1std_min = max(list_ind_less_than_1std_from_min)

    dev_allow_range = minimum_mean_dev+(0.5*std_of_min_dev)
    list_ind_less_than_05std_from_min = list(filter(lambda ind: mean_dev_per_penalizer[ind]<dev_allow_range, range(len(list_of_list_dev))))
    ind_maximum_dev_from_05std_min = max(list_ind_less_than_05std_from_min)
    
    return ind_min_mean_dev, ind_maximum_dev_from_1std_min, ind_maximum_dev_from_05std_min
