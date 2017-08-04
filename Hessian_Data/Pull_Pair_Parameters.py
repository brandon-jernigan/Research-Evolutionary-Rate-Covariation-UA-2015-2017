# -*- coding: utf-8 -*-
"""
Created on Sat Jul 16 16:48:26 2016

@author: brandon jernigan

Pull_Pair_Parameters.py this script compares a function of co-influence with the ERC values for brian's proteins

+"""

import pandas as pd
import numpy as np
import os
import pickle
import matplotlib.pyplot as plt
from itertools import combinations_with_replacement
from itertools import combinations
import time
import glob
from scipy.stats.mstats import gmean
import scipy.stats  as stats

def main():
    
    Pull_Pair_Parameters()
    
def Pull_Pair_Parameters():
    t0 = time.time()
    print("1) %ss" % (time.time() - t0))
    
#    df_erc = pd.read_csv("../ERC/yeast_tree_file_brians_domains.txt", index_col = False, na_values = [1.00, 0] , sep = '\t')
#    df_erc = df_erc.set_index(df_erc.columns)
#    print("2) %ss" % (time.time() - t0))
#    df_pval = pd.read_csv("../Streamline_ERC/streamlined_ERC_pfam/ERC/erc_pval_yeast_fixed.tsv", index_col = False, na_values = [1.00, 0] , sep = '\t')
#    df_pval = df_pval.set_index(df_pval.columns)
#    print("3) %ss" % (time.time() - t0))
    domain_locations = pickle.load( open( "domain_locations.p", "rb" ) )
    names = glob.glob(os.path.join("Systems_parameters", "*"))
    systems = []
    yeast_systems = []
    
    for name in names:
        systems.append(name.split('_')[2])
    
    systems = list(set(systems))
    print("2) %ss" % (time.time() - t0))
    
    all_dfs = {}
    
    for system in systems:
        
        path = os.path.join("Systems_parameters", "journal_%s_domains.csv" % system)
        doms = pd.read_csv(path)
       
        # make sure the number of amino acids in Brian's data matches the number in the domain annotations from pfam
        for domain in domain_locations:
            try:
                for ii, dom in enumerate(doms.loc[:]["Gene ID"]):
                    yeast_systems.append(system)
                    if dom == domain[0][0]:
    
                        for domain_num in domain[1].items():
                            match_domain = True
                            match_num = 0
                            for rng_1 in domain_num[1]:
                                
                                for rng_2 in doms.loc[ii]["Amino acids"].split():
                                    
                                    if (int(rng_1.split("-")[0]) - 1) == int(rng_2.split("-")[0]) and int(rng_1.split("-")[1]) == int(rng_2.split("-")[1].replace(",", "")):
                                        match_num += 1
                                        continue
                                    else:
                                        match_domain = False
                                        continue
                            if match_num == len(domain_num[1]):
                                match_domain = True
                            else:
                                match_domain = False
                                
                            if match_domain:
                                doms.set_value(ii, "Protein", doms.loc[ii]["Protein"].upper() + "_%s" % domain_num[0])
            except KeyError:
                pass
        doms = doms.set_index("Protein")
        all_dfs[system] = doms
        
    yeast_dfs = {}
    yeast_systems = list(set(yeast_systems))
    print("3) %ss" % (time.time() - t0))
    for yeast_system in yeast_systems:
        yeast_dfs[yeast_system] = all_dfs[yeast_system]
    
    
    domain_df_unique_list = {}
    domain_df_overlap_list = {}
    
    
    #find the reactions that both domains in the pair take part in
    for key, yeast_df in yeast_dfs.items():
        
        index = yeast_df.index.values
        domain_df_unique = pd.DataFrame(columns=index, index = index)
        domain_df_overlap = pd.DataFrame(columns=index, index = index)
        domain_df_not_unique = pd.DataFrame(columns=index, index = index)
        
        for combo in combinations_with_replacement(yeast_df.index.values, 2):
            skip = False
            
            for item in combo:
                if type(yeast_df.loc[item]["Reactions"]) == float:
                    skip = True
            if skip:
                continue
            
            reactions_list = yeast_df.loc[combo[0]]["Reactions"].split(", ") + yeast_df.loc[combo[1]]["Reactions"].split(", ")
            domain_df_not_unique = domain_df_not_unique.set_value(combo[0], combo[1], reactions_list)
            
            reactions_set = list(set(reactions_list))
            domain_df_unique.loc[combo[0]][combo[1]] = reactions_set
    
            both_list_num = [0] * len(list(set(reactions_list)))
            both_list = []
            for ii, item_1 in enumerate(list(set(reactions_list))):
                for item_2 in reactions_list:
                    if item_1 == item_2:
                        both_list_num[ii] += 1
                        if both_list_num[ii] == 2:
                            both_list.append(list(set(reactions_list))[ii])
                        
            domain_df_overlap.loc[combo[0]][combo[1]] = both_list
    
        domain_df_overlap_list[key] = domain_df_overlap 
        domain_df_unique_list[key] = domain_df_unique
        
    param_df_unique_list = {}
    print("4) %ss" % (time.time() - t0))
    
    
    #Find the parameters that pertain to the reactions found in the previous steps
    for yeast_system in yeast_systems:
        
        path = os.path.join("Systems_parameters", "journal_%s_params.csv" % yeast_system)
        params = pd.read_csv(path, index_col  = "Parameter")
    
        df_unique = domain_df_unique_list[yeast_system]
            
        index = df_unique.index.values
        param_df_unique = pd.DataFrame(columns=index, index = index)
        
        for combo in combinations_with_replacement(df_unique.index.values, 2):
            
            if type(df_unique.loc[combo[0]][combo[1]]) == float:
                continue
            for reaction in df_unique.loc[combo[0]][combo[1]]:
                for ii, rxns_associated in enumerate(params.loc[:]["Reactions"]):
                    for rxn in rxns_associated.split(", "):
                        
                        if reaction == rxn:
                            if type(param_df_unique.loc[combo[0]][combo[1]]) == float:
                                param_df_unique.loc[combo[0]][combo[1]] = []
                            param_df_unique.loc[combo[0]][combo[1]].append(params.ix[ii].name)
                            
                            param_df_unique.loc[combo[0]][combo[1]] = list(set(param_df_unique.loc[combo[0]][combo[1]]))
                        
        param_df_unique_list[yeast_system] = param_df_unique
    
    all_hess_vals = []
    all_geom_abs = []
    all_arit_abs = []
    all_arit = []
    
    geom_abs_mean_df_list = {}
    arithmetic_abs_mean_df_list = {}
    arithmetic_mean_df_list = {}
    print("5) %ss" % (time.time() - t0))
    
    #Calculate the coinfluence between pairs of domains using the Hessian matrices from Brian's previous work,
    #using values between the relevant parameters found in the previous steps.
    for yeast_system in yeast_systems:
        
        path = glob.glob(os.path.join("Hessians", "%s*_var_keys.dat" % yeast_system))[0]
        with open(path, 'r') as f_in:
            index = f_in.read().splitlines() 
            
        path = glob.glob(os.path.join("Hessians", "%s*_hessian.dat" % yeast_system))[0]
        hessians = pd.read_csv(path, sep = ' ', names = index)
        hessians.index = index
        
        df_unique = param_df_unique_list[yeast_system]
        index_2 = df_unique.index.values
        
        geom_abs_mean_df = pd.DataFrame(columns = index_2, index = index_2)
        arithmetic_abs_mean_df = pd.DataFrame(columns = index_2, index = index_2)
        arithmetic_mean_df = pd.DataFrame(columns = index_2, index = index_2)
        
        
        combo_hess_vals = []
        for combo_1 in combinations_with_replacement(df_unique.index.values, 2):
            
            hess_vals = []
            
            if type(df_unique.loc[combo_1[0]][combo_1[1]]) == float:
                continue
            
            for combo_2 in combinations(df_unique.loc[combo_1[0]][combo_1[1]], 2):
                hess_vals.append(hessians.loc[combo_2[0]][combo_2[1]])
                combo_hess_vals.append(hessians.loc[combo_2[0]][combo_2[1]])
            if combo_1[0] == "CDC14_2" and combo_1[1] == "NET1_2":
                print("CDC14_2--NET1_2")
            if combo_1[1] == "CDC14_2" and combo_1[0] == "NET1_2":
                print("CDC14_2--NET1_2")
    
            geom_abs_mean_df.loc[combo_1[0]][combo_1[1]] = gmean(np.abs(hess_vals))
            arithmetic_abs_mean_df.loc[combo_1[0]][combo_1[1]] = np.mean(np.abs(hess_vals))
            arithmetic_mean_df.loc[combo_1[0]][combo_1[1]] = np.mean(hess_vals)
        
        all_hess_vals.append(combo_hess_vals)
        all_hess_vals = all_hess_vals + combo_hess_vals
        geom_abs_mean_df_list[yeast_system] = geom_abs_mean_df
        all_geom_abs.append(geom_abs_mean_df.values.flatten())
        all_geom_abs = all_geom_abs + list(geom_abs_mean_df.values.flatten())
        arithmetic_abs_mean_df_list[yeast_system] = arithmetic_abs_mean_df
        all_arit_abs.append(arithmetic_abs_mean_df.values.flatten())
        all_arit_abs = all_arit_abs + list(arithmetic_abs_mean_df.values.flatten())
        arithmetic_mean_df_list[yeast_system] = arithmetic_mean_df
        all_arit.append(arithmetic_mean_df.values.flatten())
        all_arit = all_arit + list(arithmetic_mean_df.values.flatten())
        
    print("6) %ss" % (time.time() - t0))
    #plot all values together
    plt.hist(all_hess_vals, bins = 25, range = (-3,3), label = yeast_systems)
    plt.title("all_hess_vals")
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.show()
    plt.hist(all_geom_abs, bins = 25, range = (-3,3), label = yeast_systems)
    plt.title("all_geom_abs")
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.show()
    plt.hist(all_arit_abs, bins = 25, range = (-3,3), label = yeast_systems)
    plt.title("all_arit_abs")
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.show()
    plt.hist(all_arit, bins = 25, range = (-3,3), label = yeast_systems)
    plt.title("all_arit")
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.show()
    
    df_erc = pd.read_csv("brian_erc_fixed.tsv", index_col = 0, na_values = [1.00, 0] , sep = '\t')
    
    df_pval = pd.read_csv("brian_pval_fixed.tsv", index_col = 0, na_values = [1.00, 0] , sep = '\t')     
    
    pairs_df = pd.DataFrame(columns = ["ERC", "PVal"])
    
    #find stats, plot values for each systems biology model individually
    for yeast_system in yeast_systems:
        
        column_1 = "%s_geom_abs_mean" % yeast_system
        column_2 = "%s_arit_abs_mean" % yeast_system
        column_3 = "%s_arit_mean" % yeast_system
        pairs_df[column_1] = np.nan
        pairs_df[column_2] = np.nan
        pairs_df[column_3] = np.nan
    
    print("7) %ss" % (time.time() - t0))
    
    for yeast_system in yeast_systems: 
        for combo in combinations(geom_abs_mean_df_list[yeast_system].index, 2):
            row_name = "%s--%s" % (combo[0], combo[1])
            geom_abs_mean_col = "%s_geom_abs_mean" % yeast_system
            arit_abs_mean_col = "%s_arit_abs_mean" % yeast_system
            arit_mean_col = "%s_arit_mean" % yeast_system
            if row_name in pairs_df.index:
                pairs_df.loc[row_name][geom_abs_mean_col] = geom_abs_mean_df_list[yeast_system].loc[combo[0]][combo[1]]
                pairs_df.loc[row_name][arit_abs_mean_col] = arithmetic_abs_mean_df_list[yeast_system].loc[combo[0]][combo[1]]
                pairs_df.loc[row_name][arit_mean_col] = arithmetic_mean_df_list[yeast_system].loc[combo[0]][combo[1]]
    
                continue
            
            new_row = pd.DataFrame(columns = pairs_df.columns, index = [row_name])
            
            try:
                new_row.loc[row_name]["ERC"] = df_erc.loc[combo[0]][combo[1]]
                new_row.loc[row_name]["PVal"] = df_pval.loc[combo[0]][combo[1]]
                new_row.loc[row_name][geom_abs_mean_col] = geom_abs_mean_df_list[yeast_system].loc[combo[0]][combo[1]]
                new_row.loc[row_name][arit_abs_mean_col] = arithmetic_abs_mean_df_list[yeast_system].loc[combo[0]][combo[1]]
                new_row.loc[row_name][arit_mean_col] = arithmetic_mean_df_list[yeast_system].loc[combo[0]][combo[1]]
                
                if len(pairs_df) == 0:
                    pairs_df = new_row
                else:
                    pairs_df = pairs_df.append(new_row)
                    
            except KeyError:
                pass
    
    print("8) %ss" % (time.time() - t0))
    reaction_ERC = []
    
    for yeast_system in yeast_systems: 
        geom_abs_mean_col = "%s_geom_abs_mean" % yeast_system
        arit_abs_mean_col = "%s_arit_abs_mean" % yeast_system
        arit_mean_col = "%s_arit_mean" % yeast_system
        
        selected_mean = geom_abs_mean_col
        
        pairs_df = pairs_df.astype(np.float64)
        
        corr_df = pairs_df.dropna(subset = ["ERC", selected_mean])
        reaction_ERC = reaction_ERC + corr_df.loc[:]["ERC"].tolist()
        log_corr_df = corr_df
        log_corr_df.loc[:]["ERC"] = np.log(corr_df.loc[:]["ERC"])
        log_corr_df.loc[:][selected_mean] = np.log(corr_df.loc[:][selected_mean])
        log_corr_df = log_corr_df.dropna(subset = ["ERC", selected_mean])
        
    #    corr = stats.spearmanr(corr_df.loc[:]["ERC"], corr_df.loc[:][selected_mean])
    #    log_corr = stats.spearmanr(log_corr_df.loc[:]["ERC"], log_corr_df.loc[:][selected_mean])
    #    pairs_df.plot(kind = "scatter", x = "ERC", y = selected_mean, title = "R = %f    PVal = %f" % (corr[0], corr[1]))    
    
        corr = stats.pearsonr(corr_df.loc[:]["ERC"], corr_df.loc[:][selected_mean])
        log_corr = stats.pearsonr(log_corr_df.loc[:]["ERC"], log_corr_df.loc[:][selected_mean])
        pairs_df.plot(logy = True, kind = "scatter", x = "ERC", y = selected_mean, title = "log(Y axis)   R = %f    PVal = %f" % (corr[0], corr[1]))
    
        
        plt.show()
        
    all_erc = df_erc.values.flatten()
    all_erc = all_erc[~np.isnan(all_erc)]
    plt.hist(np.abs(reaction_ERC), bins = 10, range = (0,1))
    plt.title("reaction_ERC")
    plt.show()
    
    print("9) %ss" % (time.time() - t0))
    return

if __name__ == '__main__':
    main()