#!/usr/bin/env python3.6

import pandas as pd
import numpy as np
import scipy
from scipy import stats
import sys
import copy


"""

Programs Mines All Significant Multi-Dimensional Adjusted Results
and annotates the data with meta-results and VEP variant data

REQUIRED INPUT FILES

1. Testable_VEP_Results.txt
- All the annotated data from VEP

2. sig_multi_dim_adj_results.txt
- All the significant multi-dimensional adjusted results for a specific tissue

3. Meta Input File
- Meta Data File about the chicken samples (tells which samples are HFE vs LFE)

REQUIRED INPUT SETTINGS

1. Project Name
- Tells program how to parse the meta data file to identify the correct samples

2. Tissue
- Tells program how to parse the meta data file to identify the correct samples


"""



def create_chicken_meta_data_dict(meta_data_input_file):

    """

    Create the empty meta data dictionary for storing all
    the meta data results for analysis of the data

    : Param meta_data_input_file: Meta data file being examined

    : Return empty_meta_dict: Empty meta data dictionary to be filled


    """

    # Create lists to store variables
    project_list = []
    tissue_list = []

    # Open the input file
    input_file = open(meta_data_input_file, 'r')

    # Start looping over the file
    for line in input_file:

        # Skip the header line
        if line.startswith("Old"):
            continue

        # Deal with the rest of the data
        else:

            # Split the data
            split_data = line.split("\t")

            # Get the variables of interest
            project = split_data[2]
            tissue = split_data[3]

            # Add to the lists
            project_list.append(project)
            tissue_list.append(tissue)

    # Convert lists to sets and back to lists
    project_list = list(set(project_list))
    tissue_list = list(set(tissue_list))

    #### Create Inner Tissue Dictionary ###
    # Create empty inner dict
    inner_tissue_dict = {}

    # Loop over tissues
    for tissue in tissue_list:

        # Update the inner dictionary with empty dictionaries
        inner_tissue_dict.update({tissue: {}})
    ########################################

    ######### Meta Data Dictionary #########
    # Create a dictionary
    empty_meta_dict = {}
    
    # Loop over projects and in tissues for dictionary
    for project in project_list:

        # Make a deep copy of the
        deep_copy_inner_tissue_dict = copy.deepcopy(inner_tissue_dict)

        # Update the Meta Dictionary
        empty_meta_dict.update({project: deep_copy_inner_tissue_dict})
    #########################################

    # Close the file
    input_file.close()

    return(empty_meta_dict)



def  extract_chicken_meta_data(meta_data_input_file, empty_meta_data_dict):

    """

    Extract the meta data about the chickens with sample IDs as keys
    and feed efficiency status as values

    : Param meta_dta_input_file: Meta data file about the chickens
    : Param empty_meta_data_dict: Empty meta data dictionary to be filled

    : Return chicken_meta_dict: Filled in meta data dictionary

    """

    # Open the input file
    input_file = open(meta_data_input_file, 'r')

    # Rename Dictionary to prevent confusion
    chicken_meta_dict = empty_meta_data_dict

    # Start Looping over liens of the file
    for line in input_file:

        # Remove the new line for safety
        line = line.rstrip("\n")

        # Skip the header line
        if line.startswith("Old"):
            continue

        # Deal with actual data
        else:

            # Split the data
            split_data = line.split("\t")

            # Get the variables of interest
            correct_id = split_data[1]
            project = split_data[2]
            tissue = split_data[3]
            line = split_data[4]
            feed_status = split_data[6]

            # Update the chicken meta data dictionary
            chicken_meta_dict[project][tissue].update({correct_id: {'line': line,
                                                                     'feed_status': feed_status}})
            
    # Close the file
    input_file.close()

    return (chicken_meta_dict)


def extract_variant_information(vep_annotation_results_file_name):

    """

    Extract the variant information from the testable VEP results
    file created using VEP

    : Param vep_annotation_results_file_name: Name of file being parsed

    : Return variant_info_dict: Dictionary of variant information

    """

    # Create the dictionary to store the results
    variant_info_dict = {}

    # Open the file
    input_file = open(vep_annotation_results_file_name, 'r')

    # Loop over the file
    for line in input_file:

        # Remove the new line for safety
        line = line.rstrip("\n")

        # Skip the header line
        if line.startswith("#"):
            continue

        # Deal with Rest of Data
        else:

            # Split the lie
            split_line = line.split("\t")

            # Get variables of interest
            location = split_line[1]
            consequence = split_line[2]
            impact = split_line[3]
            gene_symbol = split_line[4]

            # Add results to dictionary
            variant_info_dict.update({location: {'consequence': consequence,
                                                 'impact': impact,
                                                 'gene_symbol': gene_symbol}})

    # Close the file
    input_file.close()

    return(variant_info_dict)


def create_sig_counter_ase_dict(sample_names):

    """

    Create counter dictionary based on the total number
    of samples

    : Param sample_names: list of all the samples

    : Return counter_sig_ase_dict; Dictionary to count the number of Sig ASE samples

    """

    # Create a dictionary to store values
    counter_sig_ase_dict = {}

    # Get the total sample counter
    total_sample_count = len(sample_names)

    # Loop over the range
    for value in range(1, total_sample_count + 1):

        # Update the dictionary with values
        counter_sig_ase_dict.update({value: 0})

    return(counter_sig_ase_dict)


def create_sample_index_dict(sample_names):

    """

    Creates a sample index dictionary, scrubbing the names
    to remove extra information about tissue source, so names
    can be used to look up meta-data results

    : Param sample_names: List of sample names from file header

    : Return sample_index_dict: Dictionary of sample names (key-index, values-names)

    """

    # Create dictionary to store names
    sample_index_dict = {}

    # Start a index counter
    index_counter = 0

    # Loop over the samples
    for sample in sample_names:

        # Split the name on the underscore
        split_sample_name = sample.split("_")

        # Take the last entry which the sample ID
        sample_id = split_sample_name[-1]

        # Add to the sample index dictionary
        sample_index_dict.update({index_counter: sample_id})

        # Add to the index counter
        index_counter += 1
        
    return (sample_index_dict)


def analyze_ase_behavior(chicken_meta_dict, project, tissue,
                         multi_dim_adj_file_name):

    """

    Analyze a multi-dimensional adjusted file with significant samples for ASE

    : Param chicken_meta_dict: Dictionary of all the meta data for the samples
    : Param project: Project name used to looking up values in the chicken meta data dictionary
    : Param tissue: Tissue being examined (used for dictionary lookup)
    : Param multi_dim_adj_file_name: Name of the file being examined

    : Return: output_file_name: Name of the results file from analysis, file needed
 
    """

    print ("Analyzing Multi-Dimensional Adjusted Samples")

    # Open the input file
    input_file = open(multi_dim_adj_file_name, 'r')

   # Create an output file for results
    output_file_name = 'tallied_meta_results.txt'

    # Open the files for writing
    output_file = open(output_file_name, 'w')

    # Write Headers Testable Data
    output_file.write("Chrom\tPos\tID\tRef\tAlt\tTotal_Biallelic\tTotal_ASE\tASE_Alleles\tASE_Alleles_HFE\tASE_Alleles_LFE\t"
                      + "Sig_ASE_HFE\tTotal_Biallelic_HFE\tSig_ASE_LFE\tTotal_Biallelic_LFE\n")

    # Start looping over the file
    for line in input_file:

        # Remove the new line from file
        line = line.rstrip("\n")
        # Hidden Tab After Last Entry (Bug in Coding)
        line = line.rstrip("\t")

        # Flag the header
        if line.startswith("#CHROM"):

            # Split the line
            line_split = line.split("\t")

            # Get sample names
            sample_names = line_split[9:]

            # Create sample index dictionary
            sample_index_dict = create_sample_index_dict(sample_names)

            # Create a Counter dictionary of ASE (use total samples count)
            counter_sig_ase_dict = create_sig_counter_ase_dict(sample_names)

        else:

            # Replaces all quotes in lines (Excels corrupts txts)
            line = line.replace('"', '')
            line = line.replace("'", "")

            # Split the line
            line_split = line.split("\t")

            # Get the variables of interest
            chromo = line_split[0]
            position = line_split[1]
            rs_id = line_split[2]
            reference = line_split[3]
            alternative = line_split[4]
            
            # Get the sample data
            samples_data = line_split[9:]

            # Create Counters
            sig_ASE_HFE = 0
            total_Biallelic_HFE = 0

            sig_ASE_LFE = 0
            total_Biallelic_LFE = 0

            total_Biallelic_Samples = 0
            total_sig_ASE_Samples = 0

            # Create lists to store directionality
            # (Convert lists to sets later to remove duplicates)
            sig_ASE_HFE_dir_list = []
            sig_ASE_LFE_dir_list = []

            # Keeps track of index for looking up sample name
            sample_index_counter = 0
            
            # Start looping over sample
            for sample in samples_data:

                # Get the sample ID name from Sample Index Dictionary (cleaned up ID)
                sample_ID = sample_index_dict[sample_index_counter]

                # Split the sample results
                sample_results = sample.split(":")

                # Get variables
                verdict = sample_results[0]
                counts = sample_results[2]
                significance_verdict = sample_results[4]

                # Filter out all non-biallelic samples
                if verdict != "Biallelic":

                    # Add to the sample index counter and continue
                    sample_index_counter += 1
                    continue

                # Biallelic Samples
                else:

                    # Get Ref and Alternative Counts
                    split_counts = counts.split(",")
                    ref_count = split_counts[0]
                    alt_count = split_counts[1]

                    # Get the Feed Efficiency status
                    feed_efficiency_status = chicken_meta_dict[project][tissue][sample_ID]['feed_status']

                    ################# High Feed Efficiency (HFE) Results ############################
                    # Add to counters based on results
                    if feed_efficiency_status == 'HFE' and significance_verdict == 'Fail':

                        # Add to specific counters
                        total_Biallelic_HFE += 1
                        total_Biallelic_Samples += 1
                        
                        # Add to the sample index counter and continue
                        sample_index_counter += 1
                        continue

                    elif feed_efficiency_status == 'HFE' and significance_verdict == 'Pass':

                        # Add to specific counters
                        total_Biallelic_HFE += 1
                        total_Biallelic_Samples += 1
                        # Sig Counters
                        sig_ASE_HFE += 1
                        total_sig_ASE_Samples += 1

                        # Get Directionality
                        # If Ref is the ASE Allele
                        if float(ref_count) > float(alt_count):
                            sig_ASE_HFE_dir_list.append("Ref")
                        # Alt Allele is the ASE Allele
                        else:
                            sig_ASE_HFE_dir_list.append("Alt")
                        
                        # Add to the sample index counter and continue
                        sample_index_counter += 1
                        continue

                    ################# Low Feed Efficiency (LFE) Results ############################
                    elif feed_efficiency_status == 'LFE' and significance_verdict == 'Fail':

                        # Add to specific counters
                        total_Biallelic_LFE += 1
                        total_Biallelic_Samples += 1
                        
                        # Add to the sample index counter and continue
                        sample_index_counter += 1
                        continue

                    elif feed_efficiency_status == 'LFE' and significance_verdict == 'Pass':

                        # Add to specific counters
                        total_Biallelic_LFE += 1
                        total_Biallelic_Samples += 1
                        # Sig Counters
                        sig_ASE_LFE += 1
                        total_sig_ASE_Samples += 1
                        

                        # Get Directionality
                        # If Ref is the ASE Allele
                        if float(ref_count) > float(alt_count):
                            sig_ASE_LFE_dir_list.append("Ref")
                        # Alt Allele is the ASE Allele
                        else:
                            sig_ASE_LFE_dir_list.append("Alt")
                        
                        # Add to the sample index counter and continue
                        sample_index_counter += 1
                        continue

                    # WARNING IN CASE OPTIONS ARE WRONG
                    else:

                        print ("Combo of feed efficiency status and significance not programmed correctly")
                        print ("Killing Program for Debugging")
                        sys.exit()

            # Writing Results to File Based on Findings
            # Validating Program is working correctly (comment out when done beta testing)

            ### Get the directionality of ASE alleles ###
            # Remove Duplicates and Convert To a List
            ASE_HFE_alleles_list = (list(set(sig_ASE_HFE_dir_list)))
            ASE_LFE_alleles_list = (list(set(sig_ASE_LFE_dir_list)))

            # Flag Empty Lists for HFE Alleles and Print NaN
            if len(ASE_HFE_alleles_list) == 0:
                ASE_HFE_alleles_str = 'NaN'
            # Non-Empty Lists Convert to a String
            else:
                ASE_HFE_alleles_str = ','.join(map(str, ASE_HFE_alleles_list))

            # Flag Empty Lists for LFE Alleles and Print NaN
            if len(ASE_LFE_alleles_list) == 0:
                ASE_LFE_alleles_str = 'NaN'
            # Non-Empty Lists Convert to a String
            else:
                ASE_LFE_alleles_str = ','.join(map(str, ASE_LFE_alleles_list))

            # Examine if the HFE vs LFE Alleles Match
            ase_alleles_list = list(set(sig_ASE_HFE_dir_list + sig_ASE_LFE_dir_list))
            ase_alleles_str = ','.join(map(str, ase_alleles_list))

            # Writes Tallies to File
            output_file.write(chromo + "\t" + position + "\t" + rs_id + "\t" + reference + "\t"
                              + alternative + "\t" + str(total_Biallelic_Samples) + "\t" + str(total_sig_ASE_Samples) + "\t"
                              + ase_alleles_str + "\t" + ASE_HFE_alleles_str + "\t" + ASE_LFE_alleles_str + "\t"
                              + str(sig_ASE_HFE) + "\t" + str(total_Biallelic_HFE) + "\t"
                              + str(sig_ASE_LFE) + "\t" + str(total_Biallelic_LFE) + "\n")
                    
    # Close the files
    input_file.close()
    output_file.close()

    return (output_file_name)

            
def create_final_annotated_results_file(merged_results_file_name, variant_info_dict):

    """

    Creates the final output file from the analysis, combing in all the variant data and adjusted
    p-value results

    : Param merged_results_file_name: Name of the results file
    : Param variant_info_dict: Dictionary of all the variant info to be merged into final data
    
    : Return None:

    """

    # Open the input file
    input_file = open(merged_results_file_name, 'r')

    # Output File Name
    output_file_name = "Annot_Multi_Dim_Results.txt"

    # Open the output file
    output_file = open(output_file_name , 'w')

    # Write the header to the file
    # Write Headers Testable Data
    output_file.write("Chrom\tPos\tID\tGene_Symbol\tConsequence\tImpact\tRef\tAlt\tTotal_Biallelic\tTotal_ASE\tASE_Alleles\tASE_Alleles_HFE\tASE_Alleles_LFE\t"
                      + "Sig_ASE_HFE\tTotal_Biallelic_HFE\tSig_ASE_LFE\tTotal_Biallelic_LFE\n")

    # Create a variant Counter (Use for looking adjusted p-value)
    variant_index_counter = 0

    # Open the input file
    for line in input_file:

        # Remove the new line for safety reasons
        line = line.rstrip("\n")

        # Skip the header line
        if line.startswith("Chrom"):
            continue

        # Deal with rest of data
        else:

            # Split the line on tab
            split_line = line.split("\t")

            # Get Variables of Interest
            chromosome = split_line[0]
            position = split_line[1]
            rs_id = split_line[2]

            # Combine Chromosome with Position
            location = chromosome + ":" + position + "-" + position

            # Get all the rest of data as one list
            rest_of_data = split_line[3:]
            rest_of_data_str = "\t".join(map(str, rest_of_data))

            # Get Results from Variant Info Dictionary
            gene_symbol_var_dict = variant_info_dict[location]['gene_symbol']
            consequence_var_dict = variant_info_dict[location]['consequence']
            impact_var_dict = variant_info_dict[location]['impact']

            # Write Results to the file
            output_file.write(chromosome + "\t" + position + "\t" + rs_id + "\t" + gene_symbol_var_dict + "\t"
                              + consequence_var_dict + "\t" + impact_var_dict + "\t" + rest_of_data_str + "\n")
            
            # Add to the variant counter
            variant_index_counter += 1



    # Close the files
    input_file.close()
    output_file.close()


##################################################################################################################################
##################################################################################################################################
##################################################################################################################################


def main():

    ##################### Input Files/Settings ##################
    meta_data_input_file = r'C:\Users\mjtom\Desktop\Final_VADT_Results_For_Paper\Masked_Genome_Results\Proportional_Hypothesis_Testing\Meta_Data_Chickens_No_Duplicates.txt'

    vep_annotation_results_file_name = 'VEP_Testable_Ab_Fat.txt'

    multi_dim_adj_file_name = 'Ab_Fat_sig_multi_dim_adj_results.txt'

    project = 'Feed Efficiency'

    tissue = 'Abdominal Fat'

    ##################################################################################################################################
    ################################################## DO NOT CHANGE BELOW ###########################################################
    ##################################################################################################################################
    
    print ("Starting Proportional Hypothesis Program")

    # Create Meta Data Dict for Storing Data
    empty_meta_data_dict = create_chicken_meta_data_dict(meta_data_input_file)

    # Extract the meta data about the chickens with IDs as dictionary keys
    chicken_meta_dict = extract_chicken_meta_data(meta_data_input_file, empty_meta_data_dict)

    # Extract Variant Information
    variant_info_dict = extract_variant_information(vep_annotation_results_file_name)

    # Examine Results of Meta Data
    merged_results_file_name = analyze_ase_behavior(chicken_meta_dict, project, tissue, multi_dim_adj_file_name)

    # Merge all Variant Info and Adjusted P-values into the Results File
    create_final_annotated_results_file(merged_results_file_name, variant_info_dict)
    
    print ("Done Running Program")

                
main()























