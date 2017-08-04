Author: Brandon Jernigan
Email: BrandonJernigan@email.arizona.edu
Date: 6/26/2017

Thesis Abstract:
Evolutionary rate covariation (ERC) is a phylogenetic measure of the evolutionary relationship between pairs of proteins. As proteins evolve over time, 
their rate of evolution (dN/dS) may vary. ERC measures how closely the evolutionary rates of two proteins match over a phylogeny. Proteins known to interact 
directly or indirectly tend to have higher ERC, because they typically experience similar evolutionary pressures within each lineage. Much is known about ERC 
at the whole protein level, but little is known at the domain level. Because individual functions of a protein are often performed by distinct domains, 
a focus on the domain level is expected to provide a clearer relationship between specific functions and ERC. Here we investigate ERC within and between domain 
families. In particular, we identify domain families with high ERC and investigate potential biochemical explanations.


Files:
Senior Thesis Final.docx is my biochemistry senior thesis on this project "Evolutionary Rate Covariation of Domain Families"

Research_Notes.pdf and Research_Notes.onepkg are 2 formats of my notes throughout the project. The .onepkg is easier to navigate, but only works with One Note.

pipeline.py is the pipeline run after inparanoid which takes its results and produces a list of phylogenetic trees
with branch lengths that the script ERC_concise.R can use to calculate ERC values. This can be split into domains or
can be for the whole protein. Each funciton in the pipeline has a description within the script.

pipeline_test_script.py is a python unit test script that ensures the major functions in pipeline.py work the way they are supposed to
even when they get altered.

paml_1.sh (and 2, 3, 4) is used to run paml after alignments, tree files, and control files are created for each protein.

ref_seq_to_gene_name.txt, refseq_to_uniprot.tsv, and swiss_to_name.txt, ensembl_to_refseq.tsv are used to translate between gene naming conventions in scripts

aaml_template.ctl is a template PAML control file used in pipeline.py

vert_genome_key.txt and yeast_genome_key.txt both give additonal information removed from the genome titles to make processing easier

vert_tree_topology.txt and yeast_tree_topology are the species trees used in PAML

.p files are intermediate files from pipeline.py to save run time

Folders:

pipeline_test: folder used to store data files for using pipeline_test_script.py, which tests pipeline.py

domain_data (see 7-17-2016 in Research_Notes): used to apply Brian's domain annotations from his previous project.

domain_type_analysis (see 2-13-2017, 4-3-2017 in Research_Notes): Analysis of ERC values between domain types

ERC: Location of ERC calculating script and results for proteins and domains.

heat_map_script (see 8-8-2016 in Research_Notes): contains script used to calculate the ERC values of 25 amino acid long
windows which can then be turned into a heat map of the ERC between two amino acid sequences

Hessian_Data (see 7-17-2016 in Research_Notes)

interactions (see 2-13-2017 in Research_Notes): script used to identify those proteins known to interact and look at their ERC values
