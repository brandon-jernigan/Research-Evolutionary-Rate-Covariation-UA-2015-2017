# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 21:01:41 2016

@author: Brandon Jernigan
"""

import os
import pickle
import unittest
from pipeline import (parser_ort, format_alignments, pull_and_add_domains_1,
                     pull_and_add_domains_2, pull_and_add_domains_2_brian, make_paml_trees,
                     make_paml_ctl, check_paml_scores, yeast_trees)
        
class pipeline_test_case(unittest.TestCase):
    

#   Test that parser_ort matches small file created by hand
    def test_parser_ort(self):
        
        file_path = os.path.join("pipeline_test", "parser_ort_check.txt")
        with open(file_path, 'r') as f_in:
            parser_ort_read = f_in.readlines()
            
        parser_ort_check = []
        
        for ii, item in enumerate(parser_ort_read):
            if ii % 2 == 0:
                keep = item
            else: 
                parser_ort_check.append("%s%s" % (keep, item))
                
        self.assertCountEqual(parser_ort(["sc", "ca"], "pipeline_test", test = True), parser_ort_check)
        
##   Test that format_alignments matches small file created by hand
    def test_format_alignments(self):
        
        file_path = os.path.join("pipeline_test", "format_alignments_check.txt")
        with open(file_path, 'r') as f_in:
            format_alignments_check = f_in.readlines()
            
        format_alignments_check.pop(0)
        
        self.assertCountEqual(format_alignments("pipeline_test", test = True), format_alignments_check)

##   Test that pull_and_add_domains_1 matches known good result
    def test_pull_and_add_domains_1(self):
        file_path = os.path.join("pipeline_test", "pull_and_add_domains_1_check.p")
        pull_and_add_domains_1_check = pickle.load( open( file_path, "rb" ) )
        
        self.assertCountEqual(pull_and_add_domains_1("pipeline_test", "swisspfam_test.txt", test = True), pull_and_add_domains_1_check)

##   Test that pull_and_add_domains_2 matches small file created by hand
    def test_pull_and_add_domains_2(self):
        file_path = os.path.join("pipeline_test", "pull_and_add_domains_2_check.txt")
        with open(file_path, 'r') as f_in:
            pull_and_add_domains_2_check = f_in.readlines()
            
        self.assertCountEqual(pull_and_add_domains_2("pipeline_test", "pull_and_add_domains_1_check.p"), pull_and_add_domains_2_check)

##   Test that pull_and_add_domains_2_brian matches small file created by hand
    def test_pull_and_add_domains_2_brian(self):
        file_path = os.path.join("pipeline_test", "pull_and_add_domains_2_brian_check.txt")
        with open(file_path, 'r') as f_in:
            pull_and_add_domains_2_brian_check = f_in.readlines()
            
        self.assertCountEqual(pull_and_add_domains_2_brian("pipeline_test"), pull_and_add_domains_2_brian_check)
        
##   Test that make_paml_trees matches small file created by hand  
    def test_make_paml_trees(self):
        file_path = os.path.join("pipeline_test", "make_paml_trees_check.txt")
        with open(file_path, 'r') as f_in:
            make_paml_trees_check = f_in.readlines()
        
        self.assertCountEqual(make_paml_trees("pipeline_test", "pipeline_test", test = True), make_paml_trees_check)

##   Test that make_paml_ctl matches small file created by hand    
    def test_make_paml_ctl(self):
        file_path = os.path.join("pipeline_test", "make_paml_ctl_check.txt")
        with open(file_path, 'r') as f_in:
            make_paml_ctl_check = f_in.readlines()
            
        strip_make_paml_ctl_check = [s.rstrip() for s in make_paml_ctl_check]
        
        self.assertCountEqual(make_paml_ctl("pipeline_test", test = True), strip_make_paml_ctl_check)

##   Test that yeast_trees matches small file created by hand   
    def test_yeast_trees(self):
        file_path = os.path.join("pipeline_test", "yeast_trees_check.txt")
        with open(file_path, 'r') as f_in:
            yeast_trees_check = f_in.readlines()
        
        self.assertCountEqual(yeast_trees("pipeline_test"), yeast_trees_check)
if __name__ == '__main__':
    unittest.main()