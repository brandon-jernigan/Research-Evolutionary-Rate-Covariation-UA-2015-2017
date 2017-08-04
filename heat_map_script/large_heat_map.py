# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 22:00:07 2016

@author: brandon jernigan

large_heat_map.py creates heat maps of proteins used in brains previous work
"""

import pandas as pd
import numpy as np
import glob
import time
import os


BRIAN_GENES = ['TPS2', 'STE2', 'STE18', 'STE12', 'CLN2', 'TPS1', 'CLB5', 'CDC20', 
         'TPI1', 'GPA1', 'STE5', 'MBP1', 'CDC5', 'PDC1', 'PFK1', 'FBA1', 
         'TEM1', 'CLB2', 'SIR2', 'PGI1', 'PGM2', 'MCM1', 'RPE1', 'PGK1', 
         'FUS3', 'SOL3', 'CDH1', 'CDC14', 'LTE1', 'TDH2', 'STE20', 'ADH1', 
         'NET1', 'RKI1', 'ENO1', 'CDC6', 'STE11', 'UGP1', 'MAD2', 'TKL1', 
         'NTH1', 'SIC1', 'SWI4', 'BCK2', 'BUB2', 'STE7', 'ESP1', 'GUT2', 
         'BAR1', 'CLN3', 'GND1', 'STE4', 'PDS1', 'SST2', 'CDC28', 'ZWF1', 
         'GPM1', 'RAP1', 'FAR1', 'CDC15', 'CDC55', 'CDC19']
         
def main():
    t0 = time.time()
    
    gene_1 = "GPA1"
    gene_2 = "FUS3" # Genes 1 and 2 set the genes to make a heat map for
    

    parser_ort(gene_1, gene_2)
    
    print("\n%d seconds" % (time.time() - t0))
   

def parser_ort(gene_1, gene_2):

    if not os.path.exists("heat_maps"):
        os.makedirs("heat_maps")
        
    increment = 25
    gene_pair = (gene_1, gene_2)
    yeast_file = "sub_3.tsv"
    aln_folder = "alignments_fixed_brian"
    
    df = pd.read_csv(yeast_file, index_col = 0, na_values = [1.00, 0] , sep = '\t')
    alignments = glob.glob(os.path.join(aln_folder, "aln_*"))
    
    gene_names = set()
    for window in df.columns:
        gene_names.add(window.split('_')[0])
    
    gene_lists = {}
    for name in gene_names:
        gene_lists[name] = []
        for window in df.columns:
            if name == window.split('_')[0]:
                gene_lists[name].append(window)
    
    df_subset_1 = df.loc[gene_lists[gene_pair[0]], gene_lists[gene_pair[1]]]
    df_subset_2 = df.loc[gene_lists[gene_pair[1]], gene_lists[gene_pair[0]]]
        
    translate_ref_seq = {}
    
    # make translation dict between gene naming conventions
    with open("ref_seq_to_gene_name.txt", 'r') as f_in1:
        lines = f_in1.readlines()
        for i, line in enumerate(lines):
            if i == 0:
                continue
            elif len(line.split()) < 2:
                continue
            else:
                translate_ref_seq[line.split()[1]] = line.split()[0]
    
    trans_alns = dict()

    #translate gene names
    for alignment in alignments:
        
        name = "%s_%s" % (alignment.split(os.pathsep)[1].split("_")[1], alignment.split(os.pathsep)[1].split("_")[2])
        end = name.find('.')
        
        try:
            trans_aln = translate_ref_seq[name[:end ]]
        except KeyError:                 
            trans_aln = name
            pass
        trans_alns[trans_aln] = alignment
        
    file_name_1 = trans_alns[gene_pair[0]]
    file_name_2 = trans_alns[gene_pair[1]]
    
    file_1 = []
    file_2 = []
    cere = False

    #find the sequence for s. cerevisiae, the model species, for
    # both genes being compared
    with open(file_name_1, 'r') as f_in_1:
        one = f_in_1.readlines()
        
        for ii, line in enumerate(one):
            if ii == 2:
                file_1.append(line.rstrip("\n"))
            elif line == "cerevisiae\n":
                cere = True
            elif cere:
                file_1.append(line.rstrip("\n"))
                cere = False
                break
            
    cere = False       
    with open(file_name_2, 'r') as f_in_2:
        two = f_in_2.readlines()
    
        for ii, line in enumerate(two):
            if ii == 2:
                file_2.append(line.rstrip("\n"))
            elif line == "cerevisiae\n":
                cere = True
            elif cere:
                file_2.append(line.rstrip("\n"))
                cere = False
                break
    
    seq_1 = []
    seq_2 = []
    first = True
    
    #make the table row and column labels
    for num, char in zip(file_1[0], file_1[1]):
        if first:
            seq_1.append("[%s]%s(%s)" % (gene_pair[0], char, num))
            first = False
        else:
            seq_1.append("%s(%s)" % (char, num))
            
    first = True
    
    for num, char in zip(file_2[0], file_2[1]):
        if first:
            seq_2.append("[%s]%s(%s)" % (gene_pair[1], char, num))
            first = False
        else:
            seq_2.append("%s(%s)" % (char, num))
    
    df_seq_1 = pd.DataFrame(columns=seq_2, index = seq_1)
    df_seq_2 = pd.DataFrame(columns=seq_1, index = seq_2)
    
    gene_list_nums = {}
    
    for key in gene_lists.keys():
        gene_list_nums[key] = []
        for item in gene_lists[key]:
            gene_list_nums[key].append(int(item.split('_')[1]))
    
    for key in gene_list_nums.keys():
        gene_list_nums[key] = sorted(gene_list_nums[key])
        
    value_1 = True
    value_2 = True

    test_nan_1 = df_subset_1.values.flatten()[~np.isnan(df_subset_1.values.flatten())]
    if len(test_nan_1) == 0:
        value_1 = False   
    test_nan_2 = df_subset_2.values.flatten()[~np.isnan(df_subset_2.values.flatten())]
    if len(test_nan_2) == 0:
        value_2 = False
    

    #fill in the table with ERC values between equal sized windows of the sequences.
    #the windows are 25 amino acids long.
    start_row = increment//2
    start_col = increment//2
    if value_1:
        for place_1 in gene_list_nums[gene_pair[0]]:
            print('.', end="")
            for place_2 in gene_list_nums[gene_pair[1]]:
                sub_row = "%s_%s" % (gene_pair[0], place_1)
                sub_col = "%s_%s" % (gene_pair[1], place_2)
                seq_row = increment * (place_1 - 1) + start_row
                seq_col = increment * (place_2 - 1) + start_col
                square_val = df_subset_1.loc[sub_row, sub_col]
                for ii in range(25):
                    for jj in range(25):
                        
                        adjust_row = jj - 12 
                        adjust_col = ii - 12
                        df_seq_1.ix[(seq_row + adjust_row),(seq_col + adjust_col)] = square_val
        
        df_seq_1.to_csv(os.path.join("heat_maps", "%s-%s_1.tsv" % (gene_pair[0], gene_pair[1])), sep = '\t')
           
    
    if value_2:
        for place_1 in gene_list_nums[gene_pair[0]]:
            print('.', end="")
            for place_2 in gene_list_nums[gene_pair[1]]:
                sub_row = "%s_%s" % (gene_pair[0], place_1)
                sub_col = "%s_%s" % (gene_pair[1], place_2)
                seq_row = increment * (place_1 - 1) + start_row
                seq_col = increment * (place_2 - 1) + start_col
                square_val = df_subset_2.loc[sub_col, sub_row]
                for ii in range(25):
                    for jj in range(25):
                        
                        adjust_row = jj - 12 
                        adjust_col = ii - 12
                        df_seq_2.ix[(seq_col + adjust_col),(seq_row + adjust_row)] = square_val

        df_seq_2.to_csv(os.path.join("heat_maps", "%s-%s_2.tsv" % (gene_pair[0], gene_pair[1])), sep = '\t')                 


if __name__ == '__main__':
    main()
