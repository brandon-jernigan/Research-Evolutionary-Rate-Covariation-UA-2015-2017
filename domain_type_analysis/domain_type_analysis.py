# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 19:49:25 2016

@author: Brandon Jernigan

Functions for running analyses on the relationship between domain family and ERC
"""

from __future__ import print_function
import csv
import time
from itertools import combinations, combinations_with_replacement, product
import os

import numpy as np
import pickle
import pandas as pd
import scipy.stats  as stats
import networkx as nx
import matplotlib.pyplot as plt
from statsmodels.sandbox.stats.multicomp import multipletests

ERC_TABLE = "yeast_domain_erc_values.tsv"
PHYLOGENETIC_TREES = "yeast_tree_file.txt"
CORRECTION = 'bonferroni'

def main():
    #erc_combinations.py finds the erc values between domain families
    erc_combinations(ERC_TABLE) 

    
    #erc_same_protein.py finds the erc values between domains on the same protein, used to compare to those between proteins
    erc_same_protein(ERC_TABLE)
    
    #compare_branch_lengths.py compares branch lengths between sets of phylogenetic trees
    compare_branch_lengths(PHYLOGENETIC_TREES)
    
    #network_domain_erc.py create a weighted networkx graph of the domain families' ERC values
    network_domain_erc()
    
    #network_domain_erc_2.py creates an inversely weighted networkx graph of the domain families' ERC values
    network_domain_erc_2()
    
    return


#erc_combinations.py finds the erc values between domain families
def erc_combinations(erc_values, test = False):
    #Correction methods for multipletests function
    #`bonferroni` : one-step correction
    #`sidak` : one-step correction
    #`holm-sidak` : step down method using Sidak adjustments
    #`holm` : step-down method using Bonferroni adjustments
    #`simes-hochberg` : step-up method  (independent)
    #`hommel` : closed method based on Simes tests (non-negative)
    #`fdr_bh` : Benjamini/Hochberg  (non-negative)
    #`fdr_by` : Benjamini/Yekutieli (negative)
    #`fdr_tsbh` : two stage fdr correction (non-negative)
    #`fdr_tsbky` : two stage fdr correction (non-negative)
    
    
    
    t0 = time.time()
    print("1) %ss" % (time.time() - t0))
    
    df_e = pd.read_csv(erc_values, index_col = False, na_values = [1.00, 0], sep = '\t')
    df_e = df_e.set_index(df_e.columns)
    print("2) %ss" % (time.time() - t0))
    pfam_domain_names_untrans = pickle.load(open("pfam_domain_names.p", "rb" ))
    
                    
    translate_ref_seq = {}
    print("3) %ss" % (time.time() - t0))
    
    #translates between gene naming conventions
    with open("ref_seq_to_gene_name.txt", 'r') as f_in:
            lines = f_in.readlines()
            for i, line in enumerate(lines):
                if i == 0:
                    continue
                elif len(line.split()) < 2:
                    continue
                else:
                    translate_ref_seq[line.split()[1]] = line.split()[0]
    pfam_domain_names_trans = {}
    print("4) %ss" % (time.time() - t0))
        
    for key, value in pfam_domain_names_untrans.items():            
        try:
            pfam_domain_names_trans[translate_ref_seq[key.split('.')[0]]] = pfam_domain_names_untrans[key]
        except KeyError:
            continue
    pfam_domain_names = pfam_domain_names_trans
                    
    domains = {}
    
    #makes list of protein domains present
    for key, value in pfam_domain_names.items():
        for domain in value:
            if domain[1] not in domains:
                domains[domain[1]] = []
            domains[domain[1]].append("%s_%s" %(key, domain[0]))
            
    domain_hist = []
    greater_domain_hist = []
    greater_domains = {}
    
    #collects those domain families that have at least 5 members
    for key, value in domains.items():
        if len(value) > 5:
            greater_domains[key] = value
            greater_domain_hist.append(len(value))
        domain_hist.append(len(value))
        
    plt.hist(greater_domain_hist)
    plt.show()
    index = greater_domains.keys()
    df_dom_mean = pd.DataFrame(np.nan, index = index, columns = index)
    df_dom_sem = pd.DataFrame(np.nan, index = index, columns = index)
    df_dom_std = pd.DataFrame(np.nan, index = index, columns = index)
    df_dom_tstat = pd.DataFrame(np.nan, index = index, columns = index)
    df_dom_pval = pd.DataFrame(np.nan, index = index, columns = index)
    df_dom_pval_corrected = pd.DataFrame(np.nan, index = index, columns = index)
    size_df = pd.DataFrame(np.nan, index = index, columns = index)
    sig_size_df = pd.DataFrame(np.nan, index = index, columns = index)
    value_dict = {}

    non_overlapping_values_combos = []
    
    
    print("5) %ss" % (time.time() - t0))
    
    #collects the erc value for each combination of ERC values and calculates stats
    for ii, combo in enumerate(combinations_with_replacement(index, 2)):
        err_flag = False
        if ii% 100 == 0:
            print("%d out of %d" % (ii, 8800))
        non_overlapping_values = pd.DataFrame(np.nan, index = index, columns = index)
        non_overlapping_values_combos = []
        
        list_1 = greater_domains[combo[0]]
        list_2 = greater_domains[combo[1]]
        
        try:
            values =  df_e.loc[list_1, list_2]
        except KeyError:
            err_flag = True
            pass
        
        for combo_2 in product(list_1, list_2):
            if combo_2[0].split('_')[0] != combo_2[1].split('_')[0]:
                if not err_flag:
                    non_overlapping_values.loc[combo_2] = values.loc[combo_2]
                    
                non_overlapping_values_combos.append(combo_2)
        value_dict[combo] = non_overlapping_values.stack()
        size_df.loc[combo] = non_overlapping_values.count().sum()
        df_dom_mean.loc[combo] = non_overlapping_values.stack().mean()
        df_dom_sem.loc[combo] = stats.sem(non_overlapping_values.stack(), axis = 0)
        df_dom_std.loc[combo] = non_overlapping_values.stack().std()
        df_dom_tstat.loc[combo], df_dom_pval.loc[combo] = stats.ttest_1samp(non_overlapping_values.stack(), popmean = 0)
        
    pickle.dump(df_dom_mean, open( "df_dom_mean.p", "wb" ))
    
    df_dom_pval_list = []
    df_dom_pval_corrected_list = []
    for combo in combinations_with_replacement(index, 2):
        df_dom_pval_list.append(df_dom_pval.loc[combo])
    
    print("6) %ss" % (time.time() - t0))
    
    reject, df_dom_pval_corrected_list, alphaS , alphaB= multipletests(df_dom_pval_list, alpha=0.05, 
                                                       method=CORRECTION, returnsorted=False)
    
    for item, combo in zip(df_dom_pval_corrected_list, combinations_with_replacement(index, 2)):
        df_dom_pval_corrected.loc[combo] = item
    sample_output = []
    #plots stats for whole data set
    plt.hist(df_dom_mean.stack())
    plt.title("All ERC Values \n ERC mean: %.3E" %np.mean(np.mean(df_dom_mean.stack())))
    plt.xlabel("ERC")
    plt.ylabel("Frequency")
    plt.show()
    print(np.mean(df_dom_mean.stack()))
    sample_output.append(np.mean(df_dom_mean.stack()))
    
    plt.hist(df_dom_sem.stack())
    plt.show()
    print("\nSEM:\n")
    print(np.mean(df_dom_sem.stack()))
    sample_output.append(np.mean(df_dom_sem.stack()))
    
    plt.hist(df_dom_std.stack())
    plt.show()
    print("\nstd:\n")
    print(np.mean(df_dom_std.stack()))
    sample_output.append(np.mean(df_dom_std.stack()))
    
    plt.hist(df_dom_tstat.stack())
    plt.show()
    print("\ntstat:\n")
    print(np.mean(df_dom_tstat.stack()))
    sample_output.append(np.mean(df_dom_tstat.stack()))
    
    plt.hist(df_dom_pval.stack())
    plt.title("PVals Uncorrected\n PVal mean: %.3E" %np.mean(df_dom_pval.stack()))
    plt.xlabel("PVal")
    plt.ylabel("Frequency")
    plt.show()
    print("\npval:\n")
    print(np.mean(df_dom_pval.stack()))
    sample_output.append(np.mean(df_dom_pval.stack()))
    
    plt.hist(df_dom_pval_corrected.stack())
    plt.title("PVals Corrected with %s\n PVal mean: %.3E" %( CORRECTION, np.mean(df_dom_pval_corrected.stack())))
    plt.xlabel("PVal")
    plt.ylabel("Frequency")
    plt.show()
    print("\npval corrected:\n")
    print(np.mean(df_dom_pval_corrected.stack()))
    sample_output.append(np.mean(df_dom_pval_corrected.stack()))
    
    df_dom_mean.transpose().to_csv("domain_type_means_new.tsv", sep = '\t')
    df_dom_std.transpose().to_csv("domain_type_stds_new.tsv", sep = '\t')
    
    significant_erc = {}
    significant_erc_list = []
    significant_erc_mean = {}
    
    #plots stats for only those domain family pairs with a significant p value after
    # the bonferroni correction.
    for combo in combinations_with_replacement(index, 2):
        if df_dom_pval_corrected.loc[combo] < 0.1:
    #    if df_dom_pval_corrected.loc[combo] < 0.1 and df_dom_mean.loc[combo] >= 0.2:
            significant_erc[combo] = ["erc mean: %.4f" % df_dom_mean.loc[combo], 
                    "sem: %.4f" % df_dom_sem.loc[combo], 
                    "std dev: %.4f" % df_dom_std.loc[combo],
                    "tstat: %.4f" % df_dom_tstat.loc[combo],
                    "corrected pval: %.E" % df_dom_pval.loc[combo]]
            sig_size_df.loc[combo] = size_df.loc[combo]
            significant_erc_list.append(df_dom_mean.loc[combo])
            significant_erc_mean[combo] = df_dom_mean.loc[combo]
            
    significant_erc_sorted_keys = sorted(significant_erc_mean.items(), key=lambda x: x[1], reverse=True)
    
    plt.hist(significant_erc_list)
    plt.title("Significant ERC Values Corrected with %s \n ERC mean: %.3E" %(CORRECTION, np.mean(significant_erc_list)))
    plt.xlabel("ERC")
    plt.ylabel("Frequency")
    plt.show()
    print("\nERC:\n")
    print(np.mean(significant_erc_list))
    
    # plot highest ERC value domain family combos
    for ii, (key, item) in enumerate(significant_erc_sorted_keys):
        
            plt.hist(value_dict[key].tolist(),bins= np.arange(-1,1,.1))
            plt.title("%s -- %s\n" % (key[0], key[1]))
    
            plt.xlim(-1,1)
            plt.xlabel("ERC Value")
            plt.ylabel("Frequency")
            plt.show()
            
            if ii >= 40:
                break
            
    # looking at the difference in number of domains within families between
    # all members and just the significant ones.
    if not test:
        print("Mean Data Number:  %d" % size_df.stack().mean())
        plt.hist(size_df.stack(), bins= np.arange(0,500,50))
        plt.title("Number of Domains Within All Families")
        plt.xlim(0,500)
        plt.xlabel("Number of Domains in Each Family")
        plt.ylabel("Frequency")
        plt.show()
        
        print("Mean Sig Data Number:  %d" % sig_size_df.stack().mean())
        plt.hist(sig_size_df.stack(), bins= np.arange(0,500,50))
        plt.title("Number of Domains for Significant, High ERC Families")
        plt.xlim(0,500)
        plt.xlabel("Number of Domains in Each Family")
        plt.ylabel("Frequency")
        plt.show()
        
    print("4) %ss" % (time.time() - t0))
    
    # output for notes
    with open('dict.csv', 'w') as csv_file:
        writer = csv.writer(csv_file)
        for key, value in significant_erc_mean.items():
           writer.writerow(["%s--%s" %(key[0], key[1]), value, df_dom_pval_corrected.loc[key], value_dict[key].size])
    
#    pickle.dump(sample_output, open("erc_combinations_check.p", 'wb'))
    return sample_output


#erc_same_protein.py finds the erc values between domains on the same protein, used to compare to those between proteins
def erc_same_protein(erc_values):
    
    t0 = time.time()
    print("1) %ss" % (time.time() - t0))
    
    sample_output = []
    
    df_e = pd.read_csv(erc_values, index_col = False, na_values = [1.00, 0] , sep = '\t')
    df_e = df_e.set_index(df_e.columns)
    proteins = {}
    protein_ERC = {}
    protien_ERC_avg = {}
    protien_ERC_std = {}
    protien_ERC_sem = {}
    protien_ERC_tstat = {}
    protien_ERC_pval = {}
    prot = ''
    
    #take out the values for Brian's proteins
    for item in df_e.columns:
        if item.split('_')[1] == 'brian':
            continue
        if prot != item.split('_')[0]:
            prot = item.split('_')[0]
            proteins[prot] = []
            protein_ERC[prot] = []
            protien_ERC_avg[prot] = []
            protien_ERC_sem[prot] = []
            protien_ERC_std[prot] = []
        proteins[prot].append(item)
    
    protien_ERC_avg_list = []
    protien_ERC_sem_list = []
    protien_ERC_std_list = []
    protien_ERC_pval_list = []
    
    #match ERC values for domains on the same protein
    for key, value in proteins.items():
        for combo in combinations(value, 2):
#            print(combo)
            protein_ERC[key].append(df_e.loc[combo[1], combo[0]])
        prot_mean = np.mean(protein_ERC[key])
        prot_sem = stats.sem(protein_ERC[key], axis=None)
        prot_std = np.std(protein_ERC[key])
        prot_tstat, prot_pval = stats.ttest_1samp(protein_ERC[key], popmean = 0)
        
        protien_ERC_avg[key] = prot_mean
        protien_ERC_sem[key] = prot_sem
        protien_ERC_std[key] = prot_std
        protien_ERC_tstat[key] = prot_tstat
        protien_ERC_pval[key] = prot_pval
        
        if ~np.isnan(prot_mean):
            protien_ERC_avg_list.append(prot_mean)
            
            protien_ERC_std_list.append(prot_std)
        if ~np.isnan(prot_sem):
            protien_ERC_sem_list.append(prot_sem)
        if ~np.isnan(prot_pval):
            protien_ERC_pval_list.append(prot_pval)
    #protien_ERC_avg_list = protien_ERC_avg_list[~np.isnan(protien_ERC_avg_list)]
    
    #plot average ERC, other stats
    plt.hist(protien_ERC_avg_list)
    plt.show()
    print(np.mean(protien_ERC_avg_list))
    sample_output.append(np.mean(protien_ERC_avg_list))
    
    plt.hist(protien_ERC_sem_list)
    plt.show()
    print(np.mean(protien_ERC_sem_list))
    sample_output.append(np.mean(protien_ERC_sem_list))
    
    plt.hist(protien_ERC_pval_list)
    plt.title("pvals")
    plt.show()
    print(np.mean(protien_ERC_pval_list))
    sample_output.append(np.mean(protien_ERC_pval_list))
    
    return sample_output


#Checks if value is a float
def RepresentsFloat(s):
    try: 
        float(s)
        return True
    except ValueError:
        return False
    
    
#compare_branch_lengths.py compares branch lengths between sets of phylogenetic trees
def compare_branch_lengths(phyologenetic_trees, test = False):
    sample_output = []

        
    with open(phyologenetic_trees) as f_in:
        tree_file = f_in.readlines()
    
    branch_lengths = []
    these_branch_lengths = []
    named_mean_branch_lengths = {}
    
    #parses newick phylogenetic tree notation to collect branch lengths across all proteins
    for item in tree_file:
        dom_name = item.split('\t')[0]
        delimited = item.split('\t')[1].replace('(', '').replace(')', '').replace(';', '')
        delimited = delimited.replace(' ', '').replace(':', ',').split(',')
        for part in delimited:
            if RepresentsFloat(part):
                branch_lengths.append(float(part))
                these_branch_lengths.append(float(part))
                
        named_mean_branch_lengths[dom_name] = np.mean(these_branch_lengths)
        these_branch_lengths = []
        
    branch_mean = np.mean(branch_lengths)
    sample_output.append(branch_mean)
    
    if not test:
        plt.hist(branch_lengths, bins= np.arange(0,1.1,0.1))
        plt.title("All Branch Lengths")
        plt.xlim(0,1)
        plt.xlabel("Branch Length")
        plt.ylabel("Frequency")
        plt.figure(1,figsize=(10000,10000)) 
        plt.savefig("Branch.png", dpi=4000)
        plt.show()
    
    #loads proteins with significant, high erc
    value_dict = pickle.load(open( "value_dict.p", "rb" ))
    significant_erc_mean = pickle.load(open( "significant_erc_mean.p", "rb" ))
    df_dom_mean = pickle.load(open( "df_dom_mean.p", "rb" ))
    
    sig_domains = {}
    print("ok")
    
    #Next portion of the code matches significant, high erc proteins with their branch lengths
    for key_1, value_1 in value_dict.items():
        for combo in combinations_with_replacement(df_dom_mean.index, 2):
            if key_1 == combo:
                try:
                    value_2 = df_dom_mean.loc[combo] 
                    if value_2 >= 0.5:
                        sig_domains[key_1] = value_1
                except KeyError:
                    continue
                
    domains = set()
    for key, value in sig_domains.items():
        for label, item in value.iteritems():
            domains.update(label)
    
    sig_branch_means = []
    
    for item in domains:
        for key, value in named_mean_branch_lengths.items():
            if item == key:
                sig_branch_means.append(value)
     
    # Next the distributions of branch lengths across all proteins, 
    # and across significant, high erc proteins so their distributions 
    # and averages can be compared     

    if not test:      
        plt.hist(sig_branch_means, bins= np.arange(0,1.1,0.1))

        plt.title("Significant, High ERC Branch Lengths")
        plt.xlim(0,1)
        plt.xlabel("Branch Length")
        plt.ylabel("Frequency")
        plt.figure(1,figsize=(10000,10000)) 
        plt.savefig("High branch.png", dpi=4000)
        plt.show()
        
    
    
    
    with open('dict2.csv', 'w') as csv_file:
        writer = csv.writer(csv_file)
        for key, value in significant_erc_mean.items():
           writer.writerow([key, value])
        
    return sample_output       


#network_domain_erc.py create a weighted networkx graph of the domain families' ERC values
def network_domain_erc():
    node_size = 100
    
    significant_erc_mean = pickle.load(open( "significant_erc_mean.p", "rb" ))
    
    G = nx.Graph()
    
    #add erc values as edges between domain families, which are the nodes.
    # This is for ERC values above 0.
    for key, value in significant_erc_mean.items():
        if value > 0:
            G.add_edge(key[0], key[1], weight = value)
    
    
    print(list(G.nodes()))
    edges,weights = zip(*nx.get_edge_attributes(G,'weight').items())
    
    pos= nx.nx_pydot.graphviz_layout(G)  
    nx.draw_networkx_nodes(G,pos,node_size = node_size)
    
    H=nx.Graph(G)
    
    # replace all the edge data with a dictionary of 
    # the single graphviz keyword "len" set to the weight
    for (u,v,w) in G.edges(data=True):
        H.add_edge(u,v,{'len':str(w)})
    pos= nx.nx_pydot.graphviz_layout(H)
    nx.draw_networkx(G,pos)
    
    
    plt.savefig("weighted_graph.png", dpi=1000) # save as png
    plt.show() # display
    
    return


#network_domain_erc_2.py creates an inversely weighted networkx graph of the domain families' ERC values
def network_domain_erc_2():
    significant_erc_mean = pickle.load(open( "significant_erc_mean.p", "rb" ))
    
    #add erc values as edges between domain families, which are the nodes.
    # This is for ERC values above 0.2.
    G=nx.Graph()
    for key, value in significant_erc_mean.items():
        if value > 0.2:
            G.add_edge(key[0], key[1], weight = 1- value)
            
    
    print(list(G.nodes()))
    print(len(list(G.edges())))
    
    
    H = nx.Graph()
    
    for (u,v,w) in G.edges(data=True):
        for key, value in w.items():
            H.add_edge(u,v,{'len':str(value)})
    pos = nx.nx_pydot.graphviz_layout(H)
    nx.draw_networkx(G,pos, node_size =.6, node_color = 'k', font_size = 8, font_color = 'r', width =.6, k=0.15,iterations=20)
    plt.axis('off')
    plt.figure(1,figsize=(10000,10000)) 
    plt.savefig("opposite_weighted_graph_protein_families_with_stragglers.png", dpi=4000) # save as png
    plt.show()
    
    return

if __name__ == '__main__':
    main()