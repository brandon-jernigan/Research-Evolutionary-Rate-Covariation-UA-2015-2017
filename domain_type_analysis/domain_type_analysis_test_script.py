# -*- coding: utf-8 -*-
"""
Created on Thu Aug 3 23:05:03 2016

@author: Brandon Jernigan
"""
    
import os
import pickle
import unittest
from domain_type_analysis import erc_combinations, erc_same_protein, compare_branch_lengths

class domain_type_analysis_test_case(unittest.TestCase):
    
##   Test that erc_combinations matches known good result
    def test_erc_combinations(self):
        
        table_path = os.path.join("domain_type_analysis_test", "domain_type_analysis_test_vals.tsv")
        check_file_path = os.path.join("domain_type_analysis_test", "erc_combinations_check.p")
        
        erc_combinations_check = pickle.load(open(check_file_path, "rb"))
            
                
        self.assertCountEqual(erc_combinations(table_path, test = True) , erc_combinations_check)
      
##   Test that erc_same_protein matches known good result
    def test_erc_same_protein(self):
        
        table_path = os.path.join("domain_type_analysis_test", "domain_type_analysis_test_vals.tsv")
        check_file_path = os.path.join("domain_type_analysis_test", "erc_same_protein_check.p")
        
        erc_same_protein_check = pickle.load(open(check_file_path, "rb"))
            
                
        self.assertCountEqual(erc_same_protein(table_path) , erc_same_protein_check)        
        
##   Test that compare_branch_lengths matches known good result
    def test_compare_branch_lengths(self):
        
        trees_path = os.path.join("domain_type_analysis_test", "yeast_tree_file_test_vals.txt")
        check_file_path = os.path.join("domain_type_analysis_test", "compare_branch_lengths_check.p")
        
        compare_branch_lengths_check = pickle.load(open(check_file_path, "rb"))
            
                
        self.assertCountEqual(compare_branch_lengths(trees_path, test = True) , compare_branch_lengths_check)
#        
        
if __name__ == '__main__':
    unittest.main()