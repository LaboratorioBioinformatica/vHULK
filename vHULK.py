#!/usr/bin/env python3
# coding: utf-8
# Edited August, 5th 2022
## This is vHULK: viral Host Unveiling Kit
# Developed by Bruno Iha and Deyvid Amgarten
# Creative commons

__version__ = "2.0.0"

# Required modules
import numpy as np
import pandas as pd
from Bio import SeqIO, SearchIO
import re
import sys
from pathlib import Path
import argparse
from tensorflow.keras.models import load_model
from scipy.special import entr, logsumexp
import csv
from sklearn.preprocessing import MinMaxScaler
from joblib import load
import os
import subprocess

# Verify, download and set databases
if not os.path.isdir('database'):
    print('Downloading databases. This will only be done in the first use.')
    os.system('wget -q http://projetos.lbi.iq.usp.br/phaghost/vHULK/database_Aug_2022.tar.gz')
    print('Extracting databases files...\n')
    if subprocess.call('tar -xf database_Aug_2022.tar.gz', shell=True) == 1:
        subprocess.call('rm database_Aug_2022.tar.gz', shell=True)
        print('Error extracting databases. Please this run script again or contact the developers.')
        quit()
    print('Databases are all set!')

# Full path for the script
here = Path(__file__).resolve().parent

def parse_arguments():
    parser = argparse.ArgumentParser( prog="vHULK.py", description="Predict phage draft genomes in metagenomic bins.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Rearrange arguments to proper groups to show up better in the help
    optionalArgs = parser._action_groups.pop()

    requiredArgs = parser.add_argument_group("Required arguments")
    requiredArgs.add_argument( "-i", "--input-dir", required=True, type=lambda p: Path(p).resolve(strict=True), 
        dest="input_dir", default="", help="Path to a folder containing metagenomic bins in .fa or .fasta format")

    optionalArgs.title = "Optional arguments"
    optionalArgs.add_argument( "-o", "--output-dir", required=False, type=lambda p: Path(p).resolve(),
        dest="output_dir", default="results", help="Location to store results in. It is created if it doesn't exist")
    optionalArgs.add_argument( "-t", "--threads", dest="threads", required=False,
        default=1, type=int, help="Number of CPU threads to be used by Prokka and hmmscan")
    optionalArgs.add_argument( "-m", "--models-dir", required=False, dest="models_dir", 
        type=lambda p: Path(p).resolve(), default=here / Path("database"), help="Path to directory where all models are stored.")
    optionalArgs.add_argument( "-p", "--prune", required=False, dest="prune",
        default=True, help="Set True to prune bad results or False to show them all")
    optionalArgs.add_argument( "--all", required=False, dest="write_all",
        action="store_true", help="Write predictions for all input bins/genomes, even if they "
        "were skipped (size filtered or hmmscan failed)")
    optionalArgs.add_argument( "-v", "--version", action="version", version=__version__)
    parser._action_groups.append(optionalArgs)
    return parser.parse_args()

#check if all database files are present
def check_databases_integrity( models_dir):

    if not models_dir.is_dir():
        print("ERROR :Missing databases directory {}\n".format(models_dir))
        print( "Please run this script again to make sure you have all input data required for v.HULK" )
        sys.exit(1)

    files_ok = True
    for f in ['model_genus', 'model_species', 
              'all_vogs_hmm_profiles.hmm', 'VOGs_header.txt',
              'list_of_hosts_genus.txt', 'list_of_hosts_species.txt', 
              'scaler_entropy_genus.joblib', 'scaler_entropy_species.joblib']:
        model_path = models_dir / Path(f)
        if not model_path.exists():
            files_ok = False
            print("ERROR: Missing {} from {}".format(model_path, models_dir))

    return files_ok

"""
Strip relevant pre and suffixes from a given fasta name
    Arguments:
        fasta_p: pathlib.Path instance: The path to the fasta file
    Return:
        bin_name: str: The name of the bin as string
"""
def get_bin_name(fasta_p):
    bin_name = fasta_p.name
    if bin_name.startswith("prokka_results_"): # Check if is prefixed with prokka_results
        bin_name = fasta_p.name.replace("prokka_results_", "")
    bin_name = bin_name.replace(fasta_p.suffix, "") # Remove .fa[a|sta] suffix
    return bin_name

"""
Run prokka for a given fasta file.
    Raises CalledProcessError if command doesn't end succesfully
    Arguments:
    fasta_in: pathlib.Path instance: Path to fasta file
    output_dir: pathlib.Path instance: Path to directory where results will
        be stored
    threads: int: Number of cpus/threads for prokka to use
"""
def run_prokka(fasta_in, output_dir, threads):
    out_prefix = get_bin_name(fasta_in)
    genome_dir = output_dir / Path(out_prefix)
    command_line = ( "prokka --kingdom Viruses --centre X --compliant --gcode 11 --cpus {} --force --quiet --prefix prokka_results_{} "
        "--fast --norrna --notrna --outdir {} --cdsrnaolap --noanno {}".format(threads, out_prefix, genome_dir, fasta_in) )
    return_code = subprocess.run(command_line, shell=True)
    return_code.check_returncode()

"""
Run hmmscan for a given fasta file
    Produces <prefix>_hmmscan.out and <prefix>_hmmscan.tbl.
    <prefix> is based on the name of the given file.
    If hmmscan runs successfully results are stored in the specified output_directory.
    Raises CalledProcessError if the command did not run successfully
    Arguments:
        fasta_in: pathlib.Path instance: Path to fasta file
        output_dir: pathlib.Path instance: Path to output directory where the two files will be written
        vogs_hmms: pathlib.Path instance: Path to all_vogs_hmm_profiles_feb2018.hmm
        threads: int: Number of threads for hmmscan
"""
def run_hmmscan(fasta_in, output_dir, threads, vog_hmms):
    bin_name = get_bin_name(fasta_in)
    path_hmm_out = output_dir / Path( "_".join([bin_name, "hmmscan.out"]) ) ## hmm stdout
    path_hmm_tblout = output_dir / Path( "_".join([bin_name, "hmmscan.tbl"]) ) ## hmm table output
    command_line = ( "hmmscan -o {} --cpu {} --tblout {} --noali {} {}".format(path_hmm_out, threads, path_hmm_tblout, vog_hmms, fasta_in))
    return_code = subprocess.run(command_line, shell=True)
    return_code.check_returncode()

"""
Transforms e-value to scale more readable to machine learning
"""
def scale_data(evalue):
    middlePoint = 9
    curveGrowth = 2
    maxValue = 10
    if evalue < 0:
        return 0
    elif evalue == 0:
        return maxValue
    evalue = (-np.log10(evalue)-middlePoint)/curveGrowth
    evalue = round( maxValue / (1 + np.exp(-evalue)), 1)
    return evalue
    
"""
Parse hmmscan tabular output to a dictionary.
    Arguments:
        hmmtable: pathlib.Path instance: Path to the hmmscan output, specified
            with hmmscan's --tblout option. Can also be str.
    Return:
        dic_genes_scores: dict: A dictionary with the gene ids as keys with
            a list of lists for all its hits. This is of the form
            { gene_id: [ [ hit id, (string), hit bias, (np.uint8)], ...], ...}
"""
def construct_gene_scores_matrix(hmmtable):
    dic_genes_scores = {}
    for gene in SearchIO.parse(hmmtable, "hmmer3-tab"):
        dic_genes_scores[gene.id] = []
        for hit in gene.hits:
            dic_genes_scores[gene.id].append( [hit.id, np.uint8(scale_data(hit.evalue))] )
    return dic_genes_scores

def predict_instance(x, model, hosts, species_or_genus, models_dir):
    pred = []
    pred.append(x) # tensorflow needs a batch of inputs to predict. 
    pred.append(x) # so we duplicate de number of inputs and then just get half of the results
    pred = np.array(model.predict(pred))[0]
    
    if not pred.any():
        name = "None"
        score = 0
    else:
        name = hosts[ np.argmax(pred) ]
        score = pred[ np.argmax(pred) ]
    
    scaler_entropy_genus = load( models_dir / Path('scaler_entropy_genus.joblib'))
    #scaler_energy_genus = load( models_dir / Path('scaler_energy_genus.joblib'))
    scaler_entropy_species = load( models_dir / Path('scaler_entropy_species.joblib'))
    #scaler_energy_species = load( models_dir / Path('scaler_energy_species.joblib'))
    
    if species_or_genus == 'genus':
        entropy = 2.6*sum(entr(pred))
        #energy = 10*scaler_energy_genus.transform( [[-logsumexp(pred)]] )[0][0]
    if species_or_genus == 'species':
        entropy = 10*scaler_entropy_species.transform( [[sum( entr(pred) )]] )[0][0]
        #energy = 10*scaler_energy_species.transform( [[-logsumexp(pred)]] )[0][0]
    
    return name, score, entropy

"""
Predict genus and species for one input
"""
def predict(x, genus_model, genus_hosts, species_model, species_hosts, prune, models_dir):
    name_genus, score_genus, entropy_genus = predict_instance(x, genus_model, genus_hosts, 'genus', models_dir)
    name_species, score_species, entropy_species = predict_instance(x, species_model, species_hosts, 'species', models_dir)

    if prune == True:
        if (score_genus < 0.3) or (entropy_genus > 5.0):
            name_genus, score_genus, entropy_genus = 'None', 0.0, 0.0
        if entropy_species > 6.3 or name_species.split('_')[0] != name_genus:
            name_species, score_species, entropy_species = 'None', 0.0, 0.0

    prediction = { "pred_genus": name_genus, 
                   "score_genus": score_genus, 
                   "entropy_genus": entropy_genus, 
                   "pred_species": name_species, 
                   "score_species": score_species, 
                   "entropy_species": entropy_species}
        
    return prediction

##########
## MAIN ##
##########

def main():
    args = parse_arguments()

    print("\n**Welcome to vHULK, a toolkit for phage host prediction!")

    input_dir = args.input_dir
    output_dir = args.output_dir
    models_dir = args.models_dir
    threads = args.threads
    prune = args.prune
    
    ################################
    ## Verify databases integrity ##
    ################################

    if check_databases_integrity( models_dir ):
        model_genus = load_model( models_dir / Path("model_genus") )
        model_species = load_model( models_dir / Path("model_species") )
        with open( models_dir / Path("list_of_hosts_genus.txt") , "r") as fin:
            genus_hosts = [line.rstrip() for line in fin]
        with open( models_dir / Path("list_of_hosts_species.txt") , "r") as fin:
            species_hosts = [line.rstrip() for line in fin]

    else:
        print( "Please run this script again to make sure you have all data required for v.HULK" )
        sys.exit(1)

    ######################################
    ## Checking input and output folder ##
    ######################################
    
    # Get a list of all input fastas from the dir
    list_bins = []
    for entry in input_dir.iterdir():
        if entry.is_file() and entry.suffix.lower() in [".fa", ".fasta"]:
            list_bins.append(entry)

    # Check if input exists
    if len(list_bins) == 0:
        print("**ERROR: Input folder appears to be empty or does not contain expected bin files (suffixed with 'fa' or '.fasta').")
        sys.exit(1)
    else:
        print("**Found {} genome files in the input folder.".format(len(list_bins)))

    if not output_dir.is_dir():
        output_dir.mkdir(parents=True, exist_ok=True)

    ############
    ## Prokka ##
    ############
    
    count_prokka = 0
    prokka_skipped = []
    prokka_dir = output_dir / Path("prokka")

    print("**Prokka has started, this may take a while.")

    for bin_fasta in list_bins:
        len_bin = 0
        bin_name = get_bin_name(bin_fasta)
        for record in SeqIO.parse(bin_fasta, "fasta"):
            len_bin += len(record.seq)
        if len_bin < 5000:
            print("**vHULK has found a genome or bin, which is too short to code proteins (< 5000 bp).\n"
                "As CDSs are an import feature for vHULK, we will be skipping this: " + bin_fasta.name)
            prokka_skipped.append(bin_name)

        run_prokka(bin_fasta, prokka_dir, threads)

        count_prokka += 1
        if count_prokka % 10 == 0:
            print("**Done with {} genomes...".format(count_prokka))

    print("\n**PROKKA finished with no errors")
    print("**{:>{width}} : {}".format("Successfully run", count_prokka, width=20))

    if len(prokka_skipped) != 0:
        print("**{:>{width}} : {}".format("Skipped", len(prokka_skipped), width=20))

    # Filter out empty protein faas
    # This might occur when prokka doesn't call any genes
    valid_faas = {}
    skipped_faas = []
    for faa in prokka_dir.glob("**/*.faa"):
        bin_name = get_bin_name(faa)
        if faa.stat().st_size != 0:
            valid_faas[bin_name] = faa
        else:
            skipped_faas.append(bin_name)

    # Check point for further execution
    # There must be something to call hmmscan on
    if len(valid_faas) == 0:
        print("**ERROR: No valid protein fastas were found in {}".format(prokka_dir))
        sys.exit(1)

    ##################
    ## HMM SEARCHES ##
    ##################

    print("**HMM scan has started, this may take a while.")

    hmmscan_dir = output_dir / Path("hmmscan")
    if not hmmscan_dir.is_dir():
        hmmscan_dir.mkdir(parents=True, exist_ok=True)

    count_hmms = 0
    vog_profiles = models_dir / Path( "all_vogs_hmm_profiles.hmm" )
    for faa in valid_faas.values():
        run_hmmscan(faa, hmmscan_dir, threads, vog_profiles)
        count_hmms += 1
        print("**Done with {} / {} HMMs\r".format(count_hmms, len(valid_faas)),end="")
    else:
        print("\n**Done with HMMscan!")

    # Initialize a dictionary that holds some info
    dic_vogs_headers = {}
    headers_fp = models_dir / Path("VOGs_header.txt")
    with open(headers_fp, "r") as fin:
        for line in fin:
            vog = line.rstrip()
            dic_vogs_headers[vog] = np.uint8(0.0)

    #########################
    ## Parsing HMM results ##
    #########################

    ######################################
    ## Making deep learning predictions ##
    ######################################
    
    print("**Starting deep learning predictions")
    
    predictions_dir = output_dir / Path("predictions")
    if not predictions_dir.is_dir():
        predictions_dir.mkdir(parents=True, exist_ok=True)
    
    for entry in valid_faas:
        hmmtable = hmmscan_dir / Path( "_".join([entry, "hmmscan.tbl"]) )
        gene_scores = construct_gene_scores_matrix(hmmtable)
        x = pd.DataFrame(index=gene_scores.keys(),columns=dic_vogs_headers.keys(),dtype=np.uint8)
        x.fillna( value=np.uint8(0.0), inplace=True )

        for gene in gene_scores:
            for each_match in gene_scores[gene]:
                x[each_match[0]][gene] = np.uint8(each_match[1])

        x = x.sum()
        x = x.values.tolist()
        x = predict(x, model_genus, genus_hosts, model_species, species_hosts, prune, models_dir)
        
        # Print results for this entry
        with open(predictions_dir / Path(entry+".csv"), 'w') as f:
            w = csv.DictWriter(f, x.keys())
            w.writeheader()
            w.writerow(x)

    if args.write_all:
        if len(prokka_skipped) != 0:
            with open(output_dir / Path("skipped_due_to_prokka.txt"), 'w') as f:
                for item in prokka_skipped:
                    f.write(f"{item}\n")
        if len(skipped_faas) != 0:
            with open(output_dir / Path("skipped_no_proteins_detected.txt"), 'w') as f:
                for item in skipped_faas:
                    f.write(f"{item}\n")

    print("**Deep learning predictions have finished. "
        "Predictions are in folder {}.\n"
        "**Thank you for using vHULK".format(output_dir))

if __name__ == "__main__":
    main()
