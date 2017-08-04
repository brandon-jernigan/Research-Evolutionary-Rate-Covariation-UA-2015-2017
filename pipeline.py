#!/usr/bin/env python

#__author__ = Brandon Jernigan
#__email__ = brandonjernigan@email.arizona.edu
#__date__ = 2016-07-04

import csv
import re
from Bio import SeqIO
from Bio import Phylo
from Bio.Seq import Seq
import pandas
import glob
import time
import pickle
import os
import shutil

SPECIES = list(["sc", "ca", "cd", "cg", "cl", "ct", "dh", "eg", "kl", "kt", "le", "mg", "ss"])
SPECIES_NAMES = list(["cerevisiae", "glabrata", "lactis", "gossypii", "albicans", "dubliniens", "tropicalis", "guilliermo", "lusitaniae", "elongispor", "hansenii", "thermotole", "stipitis"])
FOLDERS = ["one", "two", "three", "four"]
ALIGNMENT_LOCATION = "alignments"

def main():
    functions = [None] * 10
    
    #which functions are active during run
    functions[0] = False #parser_ort
    functions[1] = False #format_alignments
    functions[2] = False #pull_and_add_domains_1, pull_and_add_domains_2
    functions[3] = False #pull_and_add_domains_1, pull_and_add_domains_2_brian
    functions[4] = False #partition
    functions[5] = False #make_paml_trees
    functions[6] = False #multiply_files
    functions[7] = False #make_paml_ctl
    functions[8] = False #check_paml_scores
    functions[9] = True #yeast_trees
    

    #parser_ort parses all the orthologs from inparanoid so that a list of genes files
    #with the orthologs from the each species is produced
    
    if functions[0] == True:
        parser_ort(SPECIES,"inparanoid_yeast")
    
    #After MUSCLE has been run on the genes, files with the prefix aln_ are produced in
    #a folder called Alignments. These files are formatted in format_alignments so that 
    #PAML can be run on them later 
    if functions[1] == True:
        format_alignments(ALIGNMENT_LOCATION)

    #pull_and_add domains 1 and 2 run sequentially and are only used if we are trying to
    #look at the domain level. They take a file from PFAM with the domain annotations
    #for each gene and applies this to the aligned orthologs, producing a number code
    #at the top of the file which indicates the domain number each amino acid belongs to.
    #An alternative to pull_and_add_domains_2 is pull_and_add_domains_2_brian which uses
    #the domain annotations from the genes in Brian's previous project in dynamical influence
    
    if functions[2] == True:
        
        pull_and_add_domains_1(ALIGNMENT_LOCATION, "swisspfam")
        pull_and_add_domains_2(ALIGNMENT_LOCATION, "pfam_match.p")
        
    elif functions[3] == True:
        
        pull_and_add_domains_1(ALIGNMENT_LOCATION, "swisspfam")
        pull_and_add_domains_2_brian(ALIGNMENT_LOCATION)
        

    #partition makes 4 folders for the alignments to go into so that the job can potentially
    #be split between 4 cores
    
    if functions[4] == True:
        partition(ALIGNMENT_LOCATION)

    #make_paml_trees takes the alignments and an overall yeast tree and trims away 
    #the branches not present for each alignment (not every species is represented in every
    #set of orthoglogs)
    
    if functions[5] == True:
        make_paml_trees(ALIGNMENT_LOCATION, "trees")

    #multiply_files multiplies the trees for each alignment, and the alignments by 10
    #(it creates 10 copies of each file). This is so that PAML can be run multiple times
    #on the same information, ensuring that there isn't any variablility in the predicitons
    
    if functions[6] == True:
        multiply_files(ALIGNMENT_LOCATION)

    #make_paml_ctl makes the control file for each alignment for PAML to use, from a
    #template control file

    if functions[7] == True:
        make_paml_ctl(ALIGNMENT_LOCATION)

    #At this point run PAML on the alignments, tree files, and control files


    #Check_paml_scores ensures the likelihood score for each alignment is above 0.01
    
    if functions[8] == True:
        check_paml_scores(ALIGNMENT_LOCATION)

    #yeast_trees creates a file containing all the genes with their phylogenetic tree
    #with branch lengths that were predicted by PAML, this file then goes into the ERC script
    
    if functions[9] == True:
        yeast_trees("paml_files_within_score")




#parser_ort parses all the orthologs from inparanoid so that a list of genes files
#with the orthologs from the each species is produced
def parser_ort(species_list, ortholog_loc, test = False):
       
    
    try: 
        os.makedirs("genes")
    except OSError:
        if not os.path.isdir("genes"):
            raise
    
    t0 = time.time()

    model_species = species_list.pop(0)
    sample_output = []
    
    # go through each species pair and pull out the orthologs from each, create a
    # gene file with all species for each
    for spec in species_list:
        print(species_list)
        folder = "%s-%s" % (model_species, spec)
        model_path = os.path.join(ortholog_loc, folder, model_species)
        print(model_path)
        fasta_model = SeqIO.parse(open(model_path),'fasta')
        fasta_model_list = []
        genes = []
        
        for item in fasta_model:
            genes.append(item.id)
            fasta_model_list.append(">%s\n%s" % (model_species, item.seq))
        
        gene_num = len(genes)
        df = pandas.DataFrame(index = genes, columns = species_list)
        
        df[model_species] = fasta_model_list
        
        for jj, name in enumerate(species_list):
            print("%0.2fs Starting %s (%d of %d)" % ((time.time() - t0), name, jj, len(species_list)))
            folder = "%s-%s" % (model_species, name)
            table_name = "table.%s-%s" % (model_species, name)
            table_path = os.path.join(ortholog_loc, folder, table_name)
            print(table_path)
            ortholog_results = list(csv.reader(open(table_path, 'r'), delimiter='\t'))
            
            ortholog_results.pop(0)
            first = []
            second = []
            ortholog_pairs = []
            
            print("\t%0.2fs Finding IDs from Inparanoid results..." % (time.time() - t0))
            for line in ortholog_results:
                first = line[2].split()[0]
                second = line[3].split()[0]
                ortholog_pairs.append([first, second])
            
            for ii, gene in enumerate(df.index):
                for orth in ortholog_pairs:
                    if re.match(orth[0], gene):
                        print("\t\t%s ID match found at %d out of %d..." % (name, ii, gene_num))
                        df[name][ii] = orth[1]
                        break
                    
            print("\t%0.2fs Finding sequences from genomes..." % (time.time() - t0))
            species_path = os.path.join(ortholog_loc, folder, name)
            fasta_species = SeqIO.parse(open(species_path),'fasta')
            
            for ii, gene in enumerate(df[name]):
                if pandas.isnull(gene):
                    continue
                
                for item in fasta_species:
                    if re.match(item.id, gene):
                        print("\t\t%s Seq match found at %d out of %d..." % (name, ii, gene_num))
                        df[name][ii] = ">%s\n%s" % (name, item.seq)
                        fasta_species = SeqIO.parse(open(species_path),'fasta')
                        break
        
        print("%0.2fs Writing results to file..." % (time.time() - t0))
        
        for gene in df.iterrows():
            with open(os.path.join("genes", "%s" % gene[0]), "w") as f_out: 
                for ii in range(len(species_list)):
                    if pandas.isnull(gene[1][ii]):
                        continue
                    f_out.write("%s\n" % gene[1][ii])
                    
                    if test:
                        sample_output.append("%s\n" % gene[1][ii])
                    
        
        print("%0.2fs All done!" % (time.time() - t0))
        
    return sample_output
    
    
#After MUSCLE has been run on the genes, files with the prefix aln_ are produced in
#a folder called Alignments. These files are formatted in format_alignments so that 
#PAML can be run on them later 
def format_alignments(alignment_loc, test = False):
    t0 = time.time()

    alignments = glob.glob(os.path.join(alignment_loc, "aln_*"))
    if test:
        alignments = [os.path.join(alignment_loc, "test_aln_NP_009378.1")]
        
    
    total = len(alignments)
    sample_output = []
    # go through each alignment and correct the format
    for jj, item in enumerate(alignments):
        sequence_list = []
        species = []
        sequence = ""
        print("\t%0.2fs File %d of %d..." % ((time.time() - t0), jj + 1, total))
        with open('%s' % item, 'r') as f_in:
            lines = f_in.readlines()
            
            for ii, line in enumerate(lines):
                if ii == 0:
                    header = line            
                elif line[0].islower():
                    if sequence != "":
                        sequence_list.append(sequence)
                    cut = line.split()
                    species.append(cut.pop(0))
                    sequence = "".join(cut)
                    
                else:
                    cut = line.split()
                    sequence += "".join(cut)
                    
        sequence_list.append(sequence)
        
        if test:
            item = os.path.join(alignment_loc, "aln_NP_009378.1")
            
        with open('%s' % item, 'w') as f_out:
            f_out.write(header)
            for spec, seq in zip(species, sequence_list):
                

                f_out.write("%s\n" % spec)
                f_out.write("%s\n" % seq)
                if test:
                    sample_output.append("%s\n" % spec)
                    sample_output.append("%s\n" % seq)
                
        return sample_output
  
#pull_and_add domains 1 and 2 run sequentially and are only used if we are trying to
#look at the domain level. They take a file from PFAM with the domain annotations
#for each gene and apply this to the aligned orthologs, producing a number code
#at the top of the file which indicates the domain number each amino acid belongs to.
def pull_and_add_domains_1(alignment_loc, domain_annotations, test = False):
    
    t0 = time.time()
    aln_path = os.path.join(alignment_loc, "aln_*")
    alignments = glob.glob(aln_path)
    
    print("%ds" % (time.time() - t0))
    
    for ii, item in enumerate(alignments):
        alignments[ii] = item[len(aln_path) - 1:]
    raw_pfam = []
    section = ''
    print("%ds" % (time.time() - t0))
    
    if test:
        file_path = os.path.join(alignment_loc, domain_annotations)
    
    else:
        file_path = domain_annotations

    with open(file_path, "r") as f_in:
        for line in f_in:
            if line == "\n" and section != '':
                raw_pfam.append(section)
                section = ''
            section += line
    raw_pfam.append(section)

    
    pfam = []
    print("%ds" % (time.time() - t0))
    
    #process raw pfam file
    for item in raw_pfam:
        
        gene = []
        item = item[1:-1]
        
        for ii, line in enumerate(item.split("\n")):
            domains = []
            
            if ii == 0:
                
                gene_name = line.split()[2].split(".")[0]
                gene.append(gene_name)
            
            if ii > 0:        
                
                domain_name = line.split()[0].split(".")[0]
                range_split = line.split()
                num_dom = int(range_split[1])
                
                for ii in range(num_dom):
                    domains.append(range_split[-1 * (ii + 1)])
                
                gene.append([domain_name, domains])
                
        pfam.append(gene)
    
    raw_pfam = None
    pfam_match = []
    translate = {}
    print("%ds" % (time.time() - t0))
    
    #make translation dictionary to translate between gene naming conventions
    with open("refseq_to_uniprot.tsv", 'r') as f_in_1:
        lines = f_in_1.readlines()
        
        for ii, line in enumerate(lines):
            if ii == 0:
                continue
            elif len(line.split()) < 2:
                continue
            else:
                translate[line.split()[0]] =  line.split()[1]
    
    genes_in_dataset = []
    
    in_common = set(alignments).intersection(translate)
    
    print("%ds" % (time.time() - t0))
    
    for item in in_common:
            
            try:
                name = translate[item]
                genes_in_dataset.append([item, name])
            except KeyError:
                pass
            
    translate = None
    in_common = None
                
    print("%ds" % (time.time() - t0))     
    
    #find matches between dataset and pfam genes
    for item in pfam:
        for gene in genes_in_dataset:
            if gene[1] == item[0]:
                item[0] = gene
                pfam_match.append(item)
    pfam = None
    genes_in_dataset = None
    
    print("%ds" % (time.time() - t0))
    
    #pickel results for feed into pull_and_add_domains_2
    pickle.dump( pfam_match, open( "pfam_match.p", "wb" ) )
    print("%ds" % (time.time() - t0))
    
    sample_output = pfam_match
    
    return sample_output
    

#pull_and_add domains 1 and 2 run sequentially and are only used if we are trying to
#look at the domain level. They take a file from PFAM with the domain annotations
#for each gene and applies this to the aligned orthologs, producing a number code
#at the top of the file which indicates the domain number each amino acid belongs to.
def pull_and_add_domains_2(alignment_loc, pickle_file):
    t0 = time.time()
    
    aln_path = os.path.join(alignment_loc, "aln_*")
    alignments = glob.glob(aln_path)
    
    print("%ds" % (time.time() - t0))
    for ii, item in enumerate(alignments):
        alignments[ii] = item[len(aln_path) - 1:]
    
    print("%ds" % (time.time() - t0))
    
    pickle_path = os.path.join(alignment_loc, pickle_file)
    pfam_match = pickle.load( open( pickle_path, "rb" ) )
    
    print("%ds" % (time.time() - t0))
    new_files = []
    track_dom_types = {}
    cnt_dom = 0
    flag = 0
    cnt_seq_dom = 0
    
    # go through each protein with domain annotations, then go through each 
    # sequence amino acid by amino acid and apply domain annotations
    for item in pfam_match:        
        domains = []
        dom_cnt = 1
        for ii, entry in enumerate(item):
            if ii >= 1:
                domains.append(entry[1])
                dom_cnt +=1
        new_file = []
        with open(os.path.join(alignment_loc, "aln_%s" % item[0][0]), "r") as f_in:
            old_file = f_in.readlines()
            sc_found = False
            for line in old_file:
                if sc_found == True:
                    sc_seq = line
                    sc_found = False
                if line == "cerevisiae\n":
                    sc_found = True
            reframe = [0] * len(sc_seq)
            increment = 0
            domains_fixed = []
            for jj, char in enumerate(sc_seq):
                if char == '-':
                    increment += 1
                reframe[jj] = increment
            for domain in domains:
                sections_fixed = []
                for section in domain:
                    changing = True
                    oldest_dom_beg = int(section.split('-')[0])
                    new_dom_beg = oldest_dom_beg
                    while changing:
                        old_dom_beg = new_dom_beg
                        new_dom_beg = oldest_dom_beg + reframe[old_dom_beg - 1]
                        if reframe[old_dom_beg - 1] == reframe[new_dom_beg - 1]:
                            changing = False
            
                    changing = True
                    oldest_dom_end = int(section.split('-')[1])
                    new_dom_end = oldest_dom_end
                    while changing:
                        old_dom_end = new_dom_end
                        new_dom_end = oldest_dom_end + reframe[old_dom_end - 1]
                        if reframe[old_dom_end - 1] == reframe[new_dom_end - 1]:
                            changing = False
                    sections_fixed.append("%d-%d" %(new_dom_beg, new_dom_end))
                domains_fixed.append(sections_fixed)
            print(domains_fixed, domains, sc_seq)
                   
                                  
            seq_len = int(old_file[0].split()[1])
            new_file.append("%s G\n" % old_file[0][:-1])
            new_file.append("G %s\n" % dom_cnt)
            dom_num = 1
            dom_seq = ['1'] * seq_len
            track_dom_type = []
            for domain in domains_fixed:
                dom_num +=1
                track_dom_type.append([dom_num, item[dom_num - 1][0]])
                
            dom_num = 1
            for domain in domains_fixed:
                cnt_dom +=1
                dom_num +=1
                for section in domain:
                    
                    flag = 1
                    for jj, char in enumerate(dom_seq):
                        if jj >= int(section.split('-')[0]) - 1 and jj <= int(section.split('-')[1]) - 1:
                            dom_seq[jj] = str(dom_num)
            if flag == 1:
                flag = 0
                cnt_seq_dom +=1
            dom_seq += "\n"
            new_file.append("%s" % ''.join(dom_seq))
            for jj, line in enumerate(old_file):
                if jj >= 1:
                    new_file.append("%s" % old_file[jj])
        new_files.append(["aln_%s" % item[0][0], new_file])
        track_dom_types[item[0][0]] = track_dom_type
    if not cnt_seq_dom == 0:
        print("Number of sequences with domains: %d\nNumber of domains: %d\n Average domains per sequence with domains: %f" % (cnt_seq_dom, cnt_dom, cnt_dom/cnt_seq_dom))
    print("%ds" % (time.time() - t0))
#    
                           
                    
    if not os.path.exists("alignments_fixed"):
        os.makedirs("alignments_fixed")
        
    sample_output = []
    for item_2 in new_files:
        with open(os.path.join("alignments_fixed", item_2[0]), "w") as f_out:
            for line in item_2[1]:
                f_out.write(line)
                sample_output.append(line)
    
    print("%ds" % (time.time() - t0))

    pickle.dump(track_dom_types, open( "pfam_domain_names.p", "wb" ) )
    
    print("%ds" % (time.time() - t0))
    
    return sample_output



#An alternative to pull_and_add_domains_2 is pull_and_add_domains_2_brian which uses
#the domain annotations from the genes in Brian's previous project in dynamical influence
def pull_and_add_domains_2_brian(alignment_loc):
    t0 = time.time()
    
    aln_path = os.path.join(alignment_loc, "aln_*")
    alignments = glob.glob(aln_path)
    
    print("%ds" % (time.time() - t0))
    for ii, item in enumerate(alignments):
        alignments[ii] = item[len(aln_path) - 1:]
    
    print("%ds" % (time.time() - t0))
        

# Code for creating pickle file used to shorten execution time is commented below.
# This code extracts domain annotations from Brian's data tables

    directory = os.path.join("domain_data", "Y*")
    
    seq_files = []
    files = glob.glob('%s*' % directory)
    for item in files:
        try:    
            if item.split('.')[1] == "Mgene=3":
                seq_files.append(item)
        except:
            pass
    brian_data = []
    for item in seq_files:
        with open(item, 'r') as f_in:
            text = f_in.readlines()
    
            name = os.path.split(item)[1].split(".")[0]
    
            brian_data.append([name, text])
    
    translate = {}
    print("%ds" % (time.time() - t0))
    with open("ensembl_to_refseq.tsv", 'r') as f_in_1:
        lines = f_in_1.readlines()
        for i, line in enumerate(lines):
            if i == 0:
                continue
            elif len(line.split()) < 2:
                continue
            else:
                translate[line.split()[0]] =  line.split()[1]
    
    
    print("%ds" % (time.time() - t0))
    
    brian_data_names = []
    
    for item in brian_data:
            
            try:
                name = translate[item[0]]
                item[0] = [item[0], name]
                brian_data_names.append(item)
            except KeyError:
                pass
            
    pickle.dump( brian_data_names, open( "intermediate.p", "wb" ) )

# End of pickeled code

    brian_data_names = pickle.load( open( "intermediate.p", "rb" ) )
    domain_locations = []
    aln_files_match = []
    
    #find matches between our dataset and the proteins in brian's dataset
    for item in brian_data_names:
        directory = os.path.join(alignment_loc, "aln_%s" % item[0][1])   
        path = glob.glob('%s*' % directory)
        if path != []:
            aln_files_match.append(path)
    codon = ""
    sequence = ""
    new_file =[]
    line_2 = ''
    brian_anno_to_alns = []
    
    # go through each protein with domain annotations, then go through each 
    # sequence amino acid by amino acid and apply domain annotations
    for ff in aln_files_match:
        for item in brian_data_names:
            if os.path.split(ff[0])[1].split('.')[0] == "aln_%s" % item[0][1]:
                
                    for ii, line in enumerate(item[1]):
                        if ii >= 3 and ii % 2 == 0:
                            for jj, character in enumerate(line):
                                if character == '\n':
                                    break
                                codon += character
                                
                                if (jj + 1) % 3 ==0:
                                    if codon != "---":
                                        
                                        amino = Seq(codon).translate()
                                        codon = ''
                                    else:
                                        amino = '-'
                                        codon = ''
                                    sequence += str(amino)
                            line = sequence 
                            sequence = ''
                            line += "\n"
                            new_file.append(line)        
                                
                        else:
                            if ii == 0:
                                new_line = line.split()
                                for kk, item_2 in enumerate(new_line):
                                    if kk == 0:
                                        line_2 += "%s " % item_2
                                    elif kk == 1:
                                        line_2 += "%d " % (int(item_2)//3)##
                                        
                                    elif kk == 2:
                                        line_2 += "%s" % item_2
    
                                        
                                line = line_2
                                line += "\n"
                                line_2 = ''
                            new_file.append(line)
                            
                    sc_found = False
                    domain_sequence = ''
                    for line in new_file:
                        
                        if sc_found == True:
                            sc_found = False
                            sc_seq = line
                            for amino, dom in zip(sc_seq, item[1][2]):
                                if amino != '-':
                                    domain_sequence += dom
                        if line == "cerevisiae\n":
                            sc_found = True
    
                    domain = ''
                    compare = '1'
                    count = 0
                    domains = {}

                    for ll, character in enumerate(item[1][2]):
                        if item[0][0] == "YGL116W":
                            print("found YGL116W")
                        if compare == character:
                            continue
                        else:
                            count += 1
                            
                            if count == 1:
                                begin = ll + 1
                                compare = character
                                
                            elif count == 2 and character == '1':
                                domain = "%s-%s" %(begin, ll)
                                count = 0
                                
                                if compare not in domains:
                                    domains[compare] = []
                                    domains[compare].append(domain)                                                            
                                else:
                                    domains[compare].append(domain)   
                                    
                                compare = character
                            else:
                                domain = "%s-%s" %(begin, ll)
                                count = 1
                                begin = ll + 1
                                
                                if compare not in domains:
                                    domains[compare] = []
                                    domains[compare].append(domain)                                                            
                                else:
                                    domains[compare].append(domain) 
                                    
                                compare = character
                    domain_locations.append([item[0], domains])
                    brian_anno_to_alns.append([item[0][1], new_file])

    pickle.dump( domain_locations, open( "domain_locations.p", "wb" ) ) 
    if not os.path.exists("alignments_fixed_brian"):
        os.makedirs("alignments_fixed_brian")
        
    for item in brian_anno_to_alns:
        with open(os.path.join("alignments_fixed_brian", item[0]), "w") as f_out:
            for line in item[1]:
                f_out.write(line)
                
    sample_output = brian_anno_to_alns[0][1]
    
    print("%ds" % (time.time() - t0))
    
    return sample_output

#partition makes 4 folders for the alignments to go into so that the job can potentially
#be split between 4 cores
def partition(alignment_loc):

    glob_path = os.path.join("genes2", "aln_*")
    aln_paths = glob.glob(glob_path)
    aln_filenames = []
    for aln_path in aln_paths:
        aln_filenames.append(os.path.split(aln_path)[1])
    part = len(aln_paths) // len(FOLDERS)
    counter = part
    
    # Partition files evenly between number of specified folders
    for folder in FOLDERS:
        counter = part    
        while len(aln_filenames) != 0:
            aln_file = aln_filenames.pop(0)
            if not os.path.exists(os.path.join(alignment_loc)):
                os.makedirs(os.path.join(alignment_loc))
            if not os.path.exists(os.path.join(alignment_loc, folder)):
                os.makedirs(os.path.join(alignment_loc, folder))
            shutil.copyfile(os.path.join("genes2", aln_file), os.path.join(alignment_loc, folder, aln_file))
            counter -= 1
            if folder != "four" and counter == 0:
                break
    return
    
   
#make_paml_trees takes the alignments and an overall yeast tree and trims away 
#the branches not present for each alignment (not every species is represented in every
#set of orthoglogs)
def make_paml_trees(alignment_loc, tree_loc, test = False):
    all_species = SPECIES_NAMES
    sample_output = []
    cnt_empty = 0
    empty = []

        
    for folder in FOLDERS:

        if test:
            alignments = glob.glob(os.path.join(alignment_loc, 'aln_*'))
        else: 
            alignments = glob.glob(os.path.join(alignment_loc, folder, 'aln_*'))
            
        trees = Phylo.parse('yeast_tree_topology.txt', "newick")
        
        try: 
            os.makedirs(os.path.join("trees", folder))
        except OSError:
            if not os.path.join(tree_loc, folder):
                raise
                
        for item in alignments:
            species = []
            trees = Phylo.parse('yeast_tree_topology.txt', "newick")
            with open(item) as f_in:
                lines = f_in.readlines()
                if len(lines) == 0:
                    cnt_empty += 1
                    empty.append(item)
                    break
                    
                for line in lines:
                    if line[0].islower() and line[0] != '-':
                        species.append(line.split()[0][:10])
            
            cut_nodes = list(set(all_species).difference(species))            
            
            for tree in trees:
                for node in cut_nodes:
                    
                    tree.prune(node)
            gene_name = os.path.basename(item)[4:]
            if test:
                file_loc = os.path.join(tree_loc,"tre_%s" % gene_name)
                Phylo.write(tree, file_loc , "newick")
                
            else:
                file_loc = os.path.join(tree_loc, folder,"tre_%s" % gene_name)
                Phylo.write(tree, file_loc, "newick")
    if test:
        with open(file_loc, 'r') as f_in:
            sample_output = f_in.readlines()
    return sample_output

#multiply_files multiplies the trees for each alignment, and the alignments by 10
#(it creates 10 copies of each file). This is so that PAML can be run multiple times
#on the same information, ensuring that there isn't any variablility in the predicitons
def multiply_files(alignment_loc):
    
    for folder in FOLDERS:
        directories = (os.path.join(alignment_loc, folder, "aln_"), os.path.join("trees", folder, "tre_"))
        for directory in directories:
            
            files = glob.glob('%s*' % directory)
        
        
            for item in files:
                for ii in range(1,10):
                    shutil.copyfile(item, "%s_%s" % (item, ii))
                os.remove(item)
            
    return

#make_paml_ctl makes the control file for each alignment for PAML to use, from a
#template control file
def make_paml_ctl(alignment_loc, test = False):
    
    for folder in FOLDERS:
        
        if test:
            alignments = glob.glob(os.path.join(alignment_loc, "aln_*" ))
        else:
            alignments = glob.glob(os.path.join(alignment_loc, folder, "aln_*" ))
        
        try: 
            os.makedirs(os.path.join("control", folder))
        except OSError:
            if not os.path.isdir(os.path.join("control", folder)):
                raise
                
        with open('aaml_template.ctl', 'r') as f_in:
            control_list = f_in.readlines()
            control_options = "".join(control_list)
            
        for item in alignments:
            gene_name = os.path.basename(item)[4:]
            header = "    seqfile = alignments/%s/aln_%s * sequence data file name\n" % (folder, gene_name)
            header += "    outfile = paml_files/%s/out_%s * main result file name\n" % (folder, gene_name)
            header += "    treefile = trees/%s/tre_%s * tree structure file name\n" % (folder, gene_name)
            
            with open(os.path.join("control", folder, "ctl_%s" % gene_name), 'w') as f_out:        
                f_out.write(header + control_options)
                sample_output = header + control_options
    sample_output = sample_output.split('\n')
    sample_output = [s.rstrip() for s in sample_output]
    return sample_output
    
    
#Check_paml_scores ensures the likelihood score for each alignment is above 0.01
def check_paml_scores(alignment_loc):
    try: 
        os.makedirs(os.path.join("paml_files_not_within_score"))
    except OSError:
        if not os.path.isdir(os.path.join("paml_files_not_within_score")):
            raise
    try: 
        os.makedirs(os.path.join("paml_files_within_score"))
    except OSError:
        if not os.path.isdir(os.path.join("paml_files_within_score")):
            raise
    cnt = 0
    
    #check PAML likelihood score in each PAML output file, put passing and failing 
    #scored files in separate folders
    for folder in FOLDERS:
        PAML_files = sorted(glob.glob(os.path.join("paml_files", folder, "out_*")))
        end = len(PAML_files)//9
        for ii in range(0, end):
            if ii % 100 == 0:
                print("%d of %d complete" %(ii, end))
            scores = []
            for jj in range(1,10):
                PAML_file = PAML_files.pop(0)
                with open("%s" %PAML_file) as f_in:
                    text = f_in.readlines()
                    for line in text:
                        if line.find("lnL") != -1:
                            scores.append(float(line.split()[4]))
            for kk, score in enumerate(scores):
                cnt = 0
                score_list = scores
                score_list.pop(kk)
                for item in score_list:
                    if abs(item - score) <= 0.01:
                        cnt += 1
                    if cnt >= 2:
                        break
                if cnt >= 2:
                    break
            name = os.path.basename(PAML_file)[:-2]
            if cnt >= 2:
                shutil.copyfile(PAML_file,  os.path.join("paml_files_within_score", "%s" % name))
            else:
                shutil.copyfile(PAML_file,  os.path.join("paml_files_not_within_score", "%s" % name))
                os.remove(PAML_file)
                
    return

#yeast_trees creates a file containing all the genes with their phylogenetic tree
#with branch lengths that were predicted by PAML, this file then goes into the ERC script
def yeast_trees(paml_file_loc):
    """
    Finds PAML yeast files in PAML directory and puts them into format 
    "name, tab, Newick tree string, newline" to be read using 
    ERC_scripts_150818.R
    """


    import glob
    import os
    # make list of desired PAML files
    try: 
        os.makedirs(os.path.join("ERC"))
    except OSError:
        if not os.path.isdir(os.path.join("ERC")):
            raise
     
    PAML_files = glob.glob(os.path.join(paml_file_loc, 'out_*'))
    # open dictionary relating file names to hgnc names
    #file_to_hgnc_map = pickle.load(open('file_id_to_hgnc.pkl','r'))
    translate_ref_seq = {}
    
    with open("ref_seq_to_gene_name.txt", 'r') as f_in1:
        lines = f_in1.readlines()
        for i, line in enumerate(lines):
            if i == 0:
                continue
            elif len(line.split()) < 2:
                continue
            else:
                translate_ref_seq[line.split()[1]] = line.split()[0]
    
    
    f_out = open(os.path.join('ERC/yeast_tree_file.txt'), 'w')
    # include only yeast trees in output, decode name and re-format data
    for PAML_file in PAML_files:
    
        sample_output = []
        f_in2 = open(PAML_file,'r')
        lines = f_in2.readlines()
        #print('Print tree for file {0}'.format(PAML_file))
        tree_occurance = 0
        f_in2.close()
    
        for line in lines:
    
            if line.startswith('(') and line.endswith(';\n'):
                tree_occurance += 1
    
                if tree_occurance % 2 == 0 and len(line.split()) > 4:
                    full_name = os.path.basename(PAML_file)[4:]
                    end = full_name.find('.')
                    try:
                        name = translate_ref_seq[full_name[:end ]]
                    except KeyError:
                        f_out.write("%s_%s\t%s" % (full_name[:end ],tree_occurance//2, line))
                        continue
                        
                    if name.find('\'') != -1:
                        continue
                    f_out.write("%s_%s\t%s" % (name, tree_occurance//2, line))
                    sample_output.append("%s_%s\t%s" % (name, tree_occurance//2, line))
  
    f_out.close()
    return sample_output



if __name__ == '__main__':
    main()