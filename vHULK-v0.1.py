#!/usr/bin/python3
# coding: utf-8
# Edited May, 27th 2020
## This is vHULK: viral Host Unveiling Kit
# Developed by Deyvid Amgarten and Bruno Iha
# Creative commons

# Import required Python modules
import numpy as np
import pandas as pd
from Bio import SeqIO
import re
import sys
import os
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
import subprocess
import datetime
import argparse
import warnings
import csv
warnings.filterwarnings("ignore", category=DeprecationWarning) 
warnings.simplefilter(action='ignore', category=FutureWarning)

from time import gmtime, strftime
from tensorflow.keras.layers import Dense, Activation, LeakyReLU, ReLU
from tensorflow.keras.models import load_model
from scipy.special import entr


# Function declarations
# Run prokka
def run_prokka(binn, input_folder, threads):
	# Check the fasta format
	prefix = get_prefix(binn)
	# Filehandle where the output of prokka will be saved
	# output_prokka = open(str(prefix)+'prokka.output', mode='w')
	# Full command line for prokka
	command_line = ('prokka --kingdom Viruses --centre X --compliant --gcode 11 --cpus ' + threads + ' --force --quiet --prefix prokka_results_' + str(prefix) + ' --fast --norrna --notrna --outdir ' + input_folder + 'results/prokka/' + str(prefix) + ' --cdsrnaolap --noanno ' + input_folder + str(binn)).split()
	return_code = subprocess.call(command_line, stderr=subprocess.PIPE)
	# Check with prokka run smothly
	if return_code == 1:
		print("Prokka may not be correctly installed. Please check that.")
		quit()

# Get prefix from bins
def get_prefix(binn):
	if re.search('.fasta', binn):
		prefix = re.sub('.fasta', '', binn)
	else:
		prefix = re.sub('.fa', '', binn)
	return(prefix)


# Extract Matrix





###
### Main code
###

# Set arguments
# Modification to use argparse
parser = argparse.ArgumentParser(description='Predic phage draft genomes in metagenomic bins.')
parser.add_argument('-i', action="store", required=True, dest="input_folder", help='Path to a folder containing metagenomic bins in .fa or .fasta format (required!)')
parser.add_argument('-t', action="store", dest="threads", default='1', help='Number of CPU threads to be used by Prokka and hmmscan (default=1)')
args = parser.parse_args()


# Greeting message
print('\n**Welcome v.HULK, a toolkit for phage host prediction!\n')



# Verify databases
if not os.path.isfile('models/all_vogs_hmm_profiles_feb2018.hmm'):
	print('**Your database and models are not set. Please, run: python download_and_set_models.py \n')
	quit()


# Create Filehandle for warnings
#warnings_handle = open('marvel-warnings.txt', 'w')

# Important variables
input_folder = args.input_folder
threads = args.threads

# Fix input folder path if missing '/'
if not re.search('/$', input_folder):
	input_folder = input_folder+'/'


# Take the input folder and list all multifasta (bins) contained inside it
#print(input_folder)
list_bins_temp = os.listdir(input_folder)
list_bins = []
count_bins = 0
# Empty folder
if list_bins_temp == []:
	print('**Input folder is empty. Exiting...\n')
	quit()
else:
	for each_bin in list_bins_temp:
		if re.search('.fasta$', each_bin, re.IGNORECASE):
			list_bins.append(each_bin)
			count_bins += 1
		elif re.search('.fa$', each_bin, re.IGNORECASE):
			list_bins.append(each_bin)
			count_bins += 1

if count_bins == 0:
	print('**There is no valid genome inside the input folder (%s).\nGenome or bins should be in \'.fasta\' or \'.fa\' format.\nExiting...'%input_folder)
	quit()

print('**Arguments are OK. Checked the input folder and found %d genomes.\n' % count_bins)
print('**'+str(datetime.datetime.now()))

# Create results folder
try:
	os.stat(input_folder + 'results/')
except:
	os.mkdir(input_folder + 'results/')


#####
# PROKKA
#####
# Running prokka for all the bins multfasta files in input folder
# Perform a check in each bin, then call the execute_prokka function individually
# It may take awhile
count_prokka = 0
print('**Prokka has started, this may take awhile. Be patient.\n')
for binn in list_bins:
	# Verify bin/Genome size
	len_bin = 0
	for record in SeqIO.parse(input_folder + binn, 'fasta'):
		len_bin += len(record.seq)
	#If a bin/genome is too short, skip it
	if len_bin < 5000:
		print('**v.HULK has found a genome or bin, which is too short to code proteins (<5000pb). As CDSs are an import feature for v.HULK, we will be skipping this: '+binn)
		continue  
	run_prokka(binn, input_folder, threads)
	count_prokka += 1
	if count_prokka % 10 == 0:
		print('**Done with %d genomes...' % count_prokka)
print('**Prokka tasks have finished!\n')


####
# HMM SEARCHES
####
print('**'+str(datetime.datetime.now()))
print('**Starting HMM scan, this may take awhile. Be patient.\n')
#print(str(datetime.datetime.now()))
# Create a new results folder for hmmscan output
try:
	os.stat(input_folder + 'results/hmmscan/')
except:
	os.mkdir(input_folder + 'results/hmmscan/')

# Call HMMscan to all genomes
dic_matrices_by_genome = {}
prop_hmms_hits = {}
count_hmm = 0
for binn in list_bins:
	# Prefix for naming results
	prefix = get_prefix(binn)
	command_line_hmmscan = 'hmmscan -o ' + input_folder + 'results/hmmscan/' + prefix + '_hmmscan.out --cpu ' + threads + ' --tblout ' + input_folder + 'results/hmmscan/' + prefix + '_hmmscan.tbl --noali models/all_vogs_hmm_profiles_feb2018.hmm ' + input_folder + 'results/prokka/' + prefix + '/prokka_results_' + prefix + '.faa'
	#print(command_line_hmmscan)
    # Use -E 1 for next time running HMMscan or leave the fix down there
    # In case hmmscan returns an error - Added only because it stopped in half
	#if os.path.exists(input_folder + 'results/hmmscan/' + prefix + '_hmmscan.tbl'):
	#	continue
	try:
		subprocess.call(command_line_hmmscan, shell=True)
		# Comment line above and uncomment line below in case you want to run v.HULK without running hmmscan all over again
		#True
	except:
		print('**Error calling HMMscan:', command_line_hmmscan)
		quit()
	count_hmm += 1
	# Iteration control
	print('**Done with %d bins HMM searches...' % count_hmm)
	## Create dictionary as ref of collumns - pVOGs
	dic_vogs_headers = {}
	with open('files/VOGs_header.txt', 'r') as file2:
		for line2 in file2:
			key = re.match('(.+)\n', line2).group(1)
			dic_vogs_headers[key] = np.float32(0.0)
#
# Parse hmmscan results by gene
	num_proteins_bin = 0
	with open(input_folder + 'results/prokka/' + prefix + '/prokka_results_' + prefix + '.faa', 'r') as faa:
		for line in faa:
			if re.search('^>', line):
				num_proteins_bin += 1
				# Get gene name here
				gene_name = re.search('^>(.*)', line).group(1)
	dic_matches = {}
    
	# Parse hmmout    
	with open(input_folder + 'results/hmmscan/' + prefix + '_hmmscan.tbl', 'r') as hmmscan_out:
		dic_genes_scores = {}
		for line in hmmscan_out:
			vog=''
			gene=''
			evalue=np.float32(0.0)
			score=np.float32(0.0)
			bias=np.float32(0.0)
			if re.match('^VOG', line):
				matches = re.match('^(VOG[\d\w]+)\s+-\s+([^\s]+)[^\d]+([^\s]+)\s+([^\s]+)\s+([^\s]+)', line)
				vog = matches[1]
				gene = matches[2]
				evalue = float(matches[3])
				score = float(matches[4])
				bias = float(matches[5])
				if gene in dic_genes_scores:
					dic_genes_scores[gene].append([vog,evalue,score,bias])
				else:
					dic_genes_scores[gene] = [[vog, evalue,score,bias]]
    # Here goes the continuation
	# Create a matrix by accession
	dic_matrices_by_genome[prefix] = pd.DataFrame(index=dic_genes_scores.keys() ,columns=dic_vogs_headers.keys(), dtype=float)
	dic_matrices_by_genome[prefix].fillna(value=np.float32(0.0), inplace=True)
	# Fill in evalue values
	for gene in dic_genes_scores:
		for each_match in dic_genes_scores[gene]:
			#print(each_match[1], gene)
            # Fix for evalue values greater than 1
			if each_match[1] > 1:
				#print(each_match[1])
				each_match[1] = 1
				#print(each_match[1])
			dic_matrices_by_genome[prefix][each_match[0]][gene] = np.float32(1.0) - np.float32(each_match[1])
print('\n**HMMscan has finished.')

# Condense matrices to array by suming up columns
list_condensed_matrices = []
list_file_names = []
for matrix in dic_matrices_by_genome:
    temp = list(dic_matrices_by_genome[matrix].sum(axis = 0, skipna = True))
    list_file_names.append(matrix)
    # Parse tag
    #if re.search('^NC_.*', matrix):
    #    matrix = matrix.replace("NC_", "NC")
    #[0]accession [1]genus [2]species
    #tags = matrix.split("_")
    # For Genus
    #temp.append(tags[1])
    #temp.append(tags[0])
    # For Species
    #temp.append(tag[1]+"_"+tag[2])
    #temp.append(tag[0])
    list_condensed_matrices.append(temp)
    
# Convert to array
#import numpy as np
array = np.array(list_condensed_matrices)
#print("ARRAY-SHAPE: ", len(array))

###
# Predictions
###

print("\n**Starting deeplearning predictions...")
# load models
model_genus_relu = load_model('models/model_genus_total_fixed_relu_08mar_2020.h5', custom_objects = {"LeakyReLU":LeakyReLU, "ReLU":ReLU})
model_genus_sm = load_model('models/model_genus_total_fixed_softmax_01mar_2020.h5', custom_objects = {"LeakyReLU":LeakyReLU, "ReLU":ReLU})
model_species_relu = load_model('models/model_species_total_fixed_relu_08mar_2020.h5', custom_objects = {"LeakyReLU":LeakyReLU, "ReLU":ReLU})
model_species_sm = load_model('models/model_species_total_fixed_softmax_01mar_2020.h5', custom_objects = {"LeakyReLU":LeakyReLU, "ReLU":ReLU})

with open(input_folder+'results/results.csv', 'w') as file:
    file.write("BIN/genome,pred_genus_relu,score_genus_relu,Pred_genus_softmax,score_genus_softmax,pred_species_relu,score_species_relu,pred_species_softmax,score_species_softmax,final_prediction,entropy\n")

for i in range(0, len(array)):
    # Genus ReLu
    #print(list_file_names[i])
    pred_gen_relu = model_genus_relu.predict(np.array([array[i]]))
    #print("Genus:ReLu")
    #print(pred_gen_relu)
    position_pred_gen_relu = np.argmax(pred_gen_relu)
    if not pred_gen_relu.any():
        name_pred_gen_relu = "None"
        score_pred_gen_relu = "0"
    else:
        list_hosts_genus = [line.rstrip('\n') for line in open('files/list_hosts_genus.txt')]
        name_pred_gen_relu = list_hosts_genus[position_pred_gen_relu]
        score_pred_gen_relu = str(pred_gen_relu[0][position_pred_gen_relu])
        #print(list_hosts_genus[position_pred_gen_relu])
        #print(position_pred_gen_relu, pred_gen_relu[0][position_pred_gen_relu])

    # Genus softmax
    pred_gen_sm = model_genus_sm.predict(np.array([array[i]]))
    #print("Genus:Softmax")
    #print(pred_gen_sm)
    position_pred_gen_sm = np.argmax(pred_gen_sm)
    list_hosts_genus = [line.rstrip('\n') for line in open('files/list_hosts_genus.txt')]
    name_pred_gen_sm = list_hosts_genus[position_pred_gen_sm]
    score_pred_gen_sm = str(pred_gen_sm[0][position_pred_gen_sm])
    #print(list_hosts_genus[position_pred_gen_sm])
    #print(position_pred_gen_sm, pred_gen_sm[0][position_pred_gen_sm])

    # Species Relu
    pred_sp_relu = model_species_relu.predict(np.array([array[i]]))
    #print("Species:ReLu")
    #print(pred_sp_relu)
    position_pred_sp_relu = np.argmax(pred_sp_relu)
    if not pred_sp_relu.any():
        name_pred_sp_relu = "None"
        score_pred_sp_relu = "0"
    else:
        list_hosts_sp = [line.rstrip('\n') for line in open('files/list_hosts_species.txt')]
        #print(list_hosts_sp)
        name_pred_sp_relu = list_hosts_sp[position_pred_sp_relu]
        score_pred_sp_relu = str(pred_sp_relu[0][position_pred_sp_relu])
        #print(list_hosts_sp[position_pred_sp_relu])
        #print(position_pred_sp_relu, pred_sp_relu[0][position_pred_sp_relu])

    # Species softmax
    pred_sp_sm = model_species_sm.predict(np.array([array[i]]))
    #print("Species:Softmax")
    #print(pred_sp_sm)
    position_pred_sp_sm = np.argmax(pred_sp_sm)
    list_hosts_sp = [line.rstrip('\n') for line in open('files/list_hosts_species.txt')]
    #print(list_hosts_sp)
    name_pred_sp_sm = list_hosts_sp[position_pred_sp_sm]
    score_pred_sp_sm = str(pred_sp_sm[0][position_pred_sp_sm])
    #print(list_hosts_sp[position_pred_sp_sm])
    #print(position_pred_sp_sm, pred_sp_sm[0][position_pred_sp_sm])
    ##
    # Calculate entropy
    entropy_genus_sm =  entr(pred_gen_sm).sum(axis=1)
    #entropy_genus_sm = "{:.7f}".format(entr(pred_gen_sm).sum(axis=1))
    #
    # Apply decision tree
    #
    final_decision = "None"
    # Relu sp
    if float(score_pred_sp_relu) > 0.9:
        final_decision = name_pred_sp_relu
    # SM sp
    if float(score_pred_sp_sm) > 0.6 and name_pred_sp_sm != final_decision:
        final_decision = name_pred_sp_sm
    # Coudn't predict species
    if final_decision == 'None':
        # Put here sm sp
        if float(score_pred_sp_sm) > 0.6:
            final_decision = name_pred_sp_sm
            #relu genus
            if float(score_pred_gen_relu) >= 0.7:
                final_decision = name_pred_gen_relu
            #sm genus
            if float(score_pred_gen_sm) >= 0.5 and name_pred_gen_sm != final_decision:
                final_decision = name_pred_gen_sm
        else:
            #relu genus
            if float(score_pred_gen_relu) >= 0.9:
                final_decision = name_pred_gen_relu
            #sm genus
            if float(score_pred_gen_sm) >= 0.4 and name_pred_gen_sm != final_decision:
                final_decision = name_pred_gen_sm
    # Predicted species.
    # Verify if genus is the same
    else:
        if re.search(name_pred_gen_relu, final_decision) or re.search(name_pred_gen_sm, final_decision):
            pass
        else:
            # relu genus
            if float(score_pred_gen_relu) >= 0.9:
                final_decision = name_pred_gen_relu
            #sm genus
            if float(score_pred_gen_sm) >= 0.5 and name_pred_gen_sm != final_decision:
                final_decision = name_pred_gen_sm

    
    # Print CSV
    with open(input_folder+'results/results.csv', 'a') as file:
        file.write(list_file_names[i]+","+name_pred_gen_relu+","+score_pred_gen_relu+","+name_pred_gen_sm+","+score_pred_gen_sm+","+name_pred_sp_relu+","+score_pred_sp_relu+","+name_pred_sp_sm+","+score_pred_sp_sm+","+final_decision+","+str(entropy_genus_sm[0])+"\n")


    #print(list_file_names[i]+","+name_pred_gen_relu+":"+score_pred_gen_relu+","+name_pred_gen_sm+":"+score_pred_gen_sm+","+name_pred_sp_relu+":"+score_pred_sp_relu+","+name_pred_sp_sm+":"+score_pred_sp_sm+","+final_decision+","+str(entropy_genus_sm))
print("\n**Deep learning predictions have finished. Results are in file \"results.csv\" inside input_folder/results/.\n**Thank you for using v.HULK")
