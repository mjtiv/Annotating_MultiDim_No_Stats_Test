# Annotating_MultiDim_No_Stats_Test

Programs Mines All Significant Multi-Dimensional Adjusted Results
and annotates the data with meta-results and VEP variant data.

VEP Website: https://uswest.ensembl.org/info/docs/tools/vep/index.html

REQUIRED INPUT FILES

1. Testable_VEP_Results.txt\
-All the annotated data from VEP

2. sig_multi_dim_adj_results.txt\
-All the significant multi-dimensional adjusted results for a specific tissue

3. Meta Input File\
-Meta Data File about the chicken samples (tells which samples are HFE vs LFE)

REQUIRED INPUT SETTINGS

1. Project Name
-Tells program how to parse the meta data file to identify the correct samples

2. Tissue
-Tells program how to parse the meta data file to identify the correct samples

PROGRAM VERSIONS

VEP output can occur in two ways and so slightly different versions of code were developed
to handle this hiccup. 

Anotating_Variants_NO_TEST_V1.py - handles VEP output where RS IDs are in first column and column header is called "#Uploaded_variation"
Anotating_Variants_NO_Test_V2.py - hadles VEP output where first column is the reference allele and is called "#Allele"
