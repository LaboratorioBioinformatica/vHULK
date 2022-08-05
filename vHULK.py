#!/usr/bin/env python
# coding: utf-8
# Edited August, 5th 2022
## This is vHULK: viral Host Unveiling Kit
# Developed by Bruno Iha and Deyvid Amgarten
# Creative commons

__version__ = "2.0.0"

# Verify, download and set databases
import os
import subprocess
if not os.path.isdir('database'):
    print('Downloading databases. This will only be done in the first use.')
    os.system('wget -q http://projetos.lbi.iq.usp.br/phaghost/vHULK/database_Jun_2022.tar.xz')
    print('Extracting databases files...\n')
    if subprocess.call('tar -xf database_Jun_2022.tar.xz', shell=True) == 1:
        subprocess.call('rm database_Jun_2022.tar.xz', shell=True)
        print('Error extracting databases. Please this run script again or contact the developers.')
        quit()
    print('Databases are all set!')
    subprocess.call('rm database_Jun_2022.tar.xz', shell=True)

# Import required modules
import numpy as np
import pandas as pd
from Bio import SeqIO, SearchIO
import re
import sys
from pathlib import Path
import datetime
import argparse
from tensorflow.keras.models import load_model
from scipy.special import entr, logsumexp

# models
MODEL_FILES = [ "model_genus", "model_species" ]

# vogs hmms
VOG_PROFILES = "all_vogs_hmm_profiles_feb2018.hmm"

# host annotation
HOST_FILES = [  "list_of_hosts_genus.txt", "list_of_hosts_species.txt", "VOGs_header.txt" ]

# Full path to this script
here = Path(__file__).resolve().parent

def parse_arguments():

    parser = argparse.ArgumentParser(
        prog="vHULK.py",
        description="Predict phage draft genomes in metagenomic bins.",
        # Display default values when printing help
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Rearrange arguments to proper groups
    # to show up better in the help
    optionalArgs = parser._action_groups.pop()

    # Required arguments
    requiredArgs = parser.add_argument_group("Required arguments")
    requiredArgs.add_argument(
        "-i",
        "--input-dir",
        required=True,
        type=lambda p: Path(p).resolve(strict=True),
        dest="input_dir",
        default="",
        help="Path to a folder containing metagenomic bins in .fa or .fasta format")

    requiredArgs.add_argument( 
        "-o", 
        "--output-dir",
        required=True,
        type=lambda p: Path(p).resolve(),
        dest="output_dir",
        default="",
        help="Location to store results in. It is created if it doesn't exist")

    # Optional args
    optionalArgs.title = "Optional arguments"
    optionalArgs.add_argument(
        "-t",
        "--threads",
        dest="threads",
        required=False,
        default=1,
        type=int,
        help="Number of CPU threads to be used by Prokka and hmmscan")

    optionalArgs.add_argument(
        "-m",
        "--models-dir",
        required=False,
        dest="models_dir",
        type=lambda p: Path(p).resolve(),
        default=here / Path("database"),
        help="Path to directory where all models are stored.")

    optionalArgs.add_argument(
        "-f",
        "--files-dir",
        required=False,
        dest="files_dir",
        default=here / Path("database"),
        help="Files directory provided with vHULK")

    optionalArgs.add_argument(
        "--all",
        required=False,
        dest="write_all",
        action="store_true",
        help="Write predictions for all input bins/genomes, even if they "
        "were skipped (size filtered or hmmscan failed)")

    optionalArgs.add_argument(
        "-v",
        "--version",
        action="version",
        version=__version__)

    parser._action_groups.append(optionalArgs)

    return parser.parse_args()

"""
Strip relevant pre and suffixes from a given fasta name
    Arguments:
        fasta_p: pathlib.Path instance: The path to the fasta file
    Return:
        bin_name: str: The name of the bin as string
"""
def get_bin_name(fasta_p):
    bin_name = fasta_p.name
    fa_suffix = fasta_p.suffix
    prokka_prefix = "prokka_results_"
    # Check if it is prefixed with prokka_results
    if bin_name.startswith(prokka_prefix):
        bin_name = fasta_p.name.replace(prokka_prefix, "")
    # Remove the .fa[a|sta] suffix
    bin_name = bin_name.replace(fa_suffix, "")
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
    command_line = (
        "prokka --kingdom Viruses --centre X --compliant "
        "--gcode 11 --cpus {} --force --quiet --prefix prokka_results_{} "
        "--fast --norrna --notrna --outdir {} "
        "--cdsrnaolap --noanno {}".format(
            threads, out_prefix, genome_dir, fasta_in
        )
    )
    return_code = subprocess.run(command_line, shell=True)

    return_code.check_returncode()

"""
Run hmmscan for a given fasta file
    Produces <prefix>_hmmscan.out and <prefix>_hmmscan.tbl.
    <prefix> is based on the name of the given file.
    If hmmscan runs successfully results are stored in the specified
    output_directory.
    Raises CalledProcessError if the command did not run successfully
    Arguments:
        fasta_in: pathlib.Path instance: Path to fasta file
        output_dir: pathlib.Path instance: Path to output directory where the
            two files will be written
        vogs_hmms: pathlib.Path instance: Path to
            all_vogs_hmm_profiles_feb2018.hmm
        threads: int: Number of threads for hmmscan
    Return:
        -
"""
def run_hmmscan(fasta_in, output_dir, vogs_hmms, threads):
    bin_name = get_bin_name(fasta_in)
    # Construct the filenames and paths
    ## for hmm stdout
    base_hmm_out = "_".join([bin_name, "hmmscan.out"])
    path_hmm_out = output_dir / Path(base_hmm_out)

    ## and hmm table output
    base_hmm_tblout = "_".join([bin_name, "hmmscan.tbl"])
    path_hmm_tblout = output_dir / Path(base_hmm_tblout)

    command_line = (
        "hmmscan -o {} --cpu {} --tblout {} --noali {} "
        "{}".format(
            path_hmm_out, threads, path_hmm_tblout, vogs_hmms, fasta_in
        )
    )

    return_code = subprocess.run(command_line, shell=True)

    return_code.check_returncode()

"""
Parse hmmscan tabular output to a dictionary.
    Arguments:
        hmmtable: pathlib.Path instance: Path to the hmmscan output, specified
            with hmmscan's --tblout option. Can also be str.
    Return:
        dic_genes_scores: dict: A dictionary with the gene ids as keys with
            a list of lists for all its hits. This is of the form
            { gene_id: [
                [ hit id, (<- string)
                  hit E-value, (<- np.float32)
                  hit bit-score, (<-np.float32)
                  hit bias, (<-np.float32)
                  ], ...],
                  ...}
"""
def construct_gene_scores_matrix(hmmtable):
    dic_genes_scores = {}
    for gene in SearchIO.parse(hmmtable, "hmmer3-tab"):
        dic_genes_scores[gene.id] = []
        for hit in gene.hits:
            hit_info = [hit.id,
                        np.float32(hit.evalue),
                        np.float32(hit.bitscore),
                        np.float32(hit.bias)]
            dic_genes_scores[gene.id].append(hit_info)
    return dic_genes_scores

"""
Predict genus
"""
def predict_genus(x, model_genus, genus_hosts):
    pred_genus = np.array(model_genus(x))
    idx_pred_genus = np.argmax(pred_genus)

    if not pred_genus.any():
        name_genus = "None"
        score_genus = 0
    else:
        name_genus = genus_hosts[idx_pred_genus]
        score_genus = pred_genus[0][idx_pred_genus]

    # Calculate entropy and energy
    entropy_genus = 2.6*sum(entr(pred_genus)[0])
    threshold_energy_genus = 4.36588
    energy_genus = 1100*(-logsumexp(pred_genus)+threshold_energy_genus)

    return name_genus, score_genus, entropy_genus, energy_genus

"""
Predict species
"""
def predict_species(x, model_species, species_hosts):
    
    pred_species = np.array(model_species(x))
    idx_pred_species = np.argmax(pred_species)
    
    if not pred_species.any():
        name_species = "None"
        score_species = 0
    else:
        name_species = species_hosts[idx_pred_species]
        score_species = pred_species[0][idx_pred_species]

    # Calculate entropy and energy
    entropy_species = 27*sum(entr(pred_species))[0]
    threshold_energy_species = 5.2379565
    energy_species = 71*(-logsumexp(pred_species)+threshold_energy_species)

    return name_species, score_species, entropy_species, energy_species

"""
Predict genus and species for all input for all models
"""
def predict(scores_array, names_array, model_genus, genus_hosts, model_species, species_hosts):
    predictions = {}

    for i in range(0, len(scores_array)):
        bin_name = names_array[i]
        bin_array = np.array([scores_array[i]], dtype=np.uint8)
        ( name_genus, score_genus, entropy_genus, energy_genus) = predict_genus(bin_array, model_genus, genus_hosts)

        ( name_species, score_species, entropy_species, energy_species) = predict_species(bin_array, model_species, species_hosts)

        predictions[bin_name] = {
            "pred_genus": name_genus,
            "score_genus": score_genus,
            "entropy_genus": entropy_genus,
            "energy_genus": energy_genus,
            "pred_species": name_species,
            "score_species": score_species,
            "entropy_species": entropy_species,
            "energy_species": energy_species}

    return predictions

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

##########
## MAIN ##
##########

def main():
    args = parse_arguments()

    print("\n**Welcome to vHULK, a toolkit for phage host prediction!")

    ## Important variables
    input_dir = args.input_dir
    output_dir = args.output_dir
    models_dir = args.models_dir
    files_dir = args.files_dir
    threads = args.threads

    if not models_dir.is_dir():
        print("ERROR :Missing models directory {}\n".format(models_dir))
        print( "Please run 'python download_and_set_models.py' \nto make sure you have all input data required for v.HULK\n" )
        sys.exit(1)

    # Verify databases integrity
    files_ok = True
    for m in MODEL_FILES:
        model_path = models_dir / Path(m)
        if not model_path.exists():
            files_ok = False
            print("ERROR: Missing {} from {}".format(model_path, models_dir))

    vog_profiles = models_dir / Path(VOG_PROFILES)
    if not vog_profiles.exists():
        files_ok = False
        print("ERROR: Missing {} from {}".format(vog_profiles, models_dir))

    for h in HOST_FILES:
        txt_fp = files_dir / Path(h)
        if not txt_fp.exists():
            files_ok = False
            print("ERROR: Missing {} from {}".format(txt_fp, files_dir))

    if not files_ok:
        print( "Please run 'python download_and_set_models.py' \nto make sure you have all input data required for v.HULK\n" )
        sys.exit(1)

    # Get a list of all input fastas from the dir
    list_bins = []
    for entry in input_dir.iterdir():
        if entry.is_file() and entry.suffix.lower() in [".fa", ".fasta"]:
            list_bins.append(entry)

    # Check if input exists
    if len(list_bins) == 0:
        print("**ERROR: Input folder appears to be empty or does not"
            "contain expected bin files (suffixed with 'fa' or '.fasta').")
        sys.exit(1)
    else:
        print("**Checked the input folder and found {} genome files.".format(len(list_bins)))

    # Check output dir and create it if needed
    if not output_dir.is_dir():
        output_dir.mkdir(parents=True, exist_ok=True)

    ############
    ## PROKKA ##
    ############
    
    count_prokka = 0
    prokka_skipped = {}
    prokka_dir = output_dir / Path("prokka")

    print("**Prokka has started, this may take a while. Be patient.")

    for bin_fasta in list_bins:
        len_bin = 0
        bin_name = get_bin_name(bin_fasta)
        for record in SeqIO.parse(bin_fasta, "fasta"):
            len_bin += len(record.seq)
        if len_bin < 5000:
            print("**vHULK has found a genome or bin, which is too short to "
                "code proteins (< 5000 bp). As CDSs are an import feature for "
                "vHULK, we will be skipping this: " + bin_fasta.name)
            prokka_skipped[bin_name] = bin_fasta
            continue

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
    # TO DO
    # Check with prokka if this is desired behavior
    valid_faas, skipped_faas = {}, {}
    for faa in prokka_dir.glob("**/*.faa"):
        bin_name = get_bin_name(faa)
        if faa.stat().st_size != 0:
            valid_faas[bin_name] = faa
        else:
            skipped_faas[bin_name] = faa

    # Check point for further execution
    # There must be something to call hmmscan on
    if len(valid_faas) == 0:
        print("**ERROR: No valid protein fastas were found in {}".format(prokka_dir))
        sys.exit(1)

    ##################
    ## HMM SEARCHES ##
    ##################

    print("**Starting HMM scan, this may take a while...")

    hmmscan_dir = output_dir / Path("hmmscan")
    if not hmmscan_dir.is_dir():
        hmmscan_dir.mkdir(parents=True, exist_ok=True)

    count_hmms = 0
    for faa in valid_faas.values():
        run_hmmscan(faa, hmmscan_dir, vog_profiles, threads)
        count_hmms += 1
        print("**Done with {} / {} HMMs\r".format(count_hmms, len(valid_faas)),end="")
    else:
        print("\n**Done with HMMscan!")

    # Initialize a dictionary that holds some info
    dic_vogs_headers = {}
    headers_fp = files_dir / Path("VOGs_header.txt")
    with open(headers_fp, "r") as fin:
        for line in fin:
            vog = line.rstrip()
            dic_vogs_headers[vog] = np.uint8(0.0)

    print("**Parsing HMM results")
    # PARSE HMM RESULTS
    dic_matrices_by_genome = {}
    for entry in valid_faas:
        hmmtbl_str = "_".join([entry, "hmmscan.tbl"])
        hmmtable = hmmscan_dir / Path(hmmtbl_str)
        gene_scores = construct_gene_scores_matrix(hmmtable)
        dic_matrices_by_genome[entry] = pd.DataFrame(index=gene_scores.keys(),columns=dic_vogs_headers.keys(),dtype=float)

        dic_matrices_by_genome[entry].fillna( value=np.uint8(0.0), inplace=True )
        for gene in gene_scores:
            for each_match in gene_scores[gene]:
                if each_match[1] > 1:
                    each_match[1] = 1
                dic_matrices_by_genome[entry][each_match[0]][gene] = np.uint8(scale_data(each_match[1]))

    list_condensed_matrices = []
    bin_names = []
    for matrix in dic_matrices_by_genome:
        temp = list(dic_matrices_by_genome[matrix].sum(axis=0, skipna=True))
        bin_names.append(matrix)
        list_condensed_matrices.append(temp)

    mat_array = np.array(list_condensed_matrices,dtype=np.uint8)
    
    #################
    ## PREDICTIONS ##
    #################
    print("**Starting deep learning predictions")

    # genus predictor
    model_genus = models_dir / Path("model_genus")
    model_genus = load_model(str(model_genus))

    # species predictor
    model_species = models_dir / Path("model_species")
    model_species = load_model(str(model_species))
    
    ## Genus hosts
    genus_hosts_fp = files_dir / Path("list_of_hosts_genus.txt")
    with open(genus_hosts_fp, "r") as fin:
        genus_hosts = [line.rstrip() for line in fin]

    ## Species hosts
    species_hosts_fp = files_dir / Path("list_of_hosts_species.txt")
    with open(species_hosts_fp, "r") as fin:
        species_hosts = [line.rstrip() for line in fin]

    # Predict and store in a dictionary
    results = predict(mat_array, bin_names, model_genus, genus_hosts, model_species, species_hosts)

    # Put it in a dataframe and write it
    results = pd.DataFrame.from_dict(results, orient="index")

    if args.write_all:
        if len(prokka_skipped) != 0 or len(skipped_faas) != 0:
            bins_skipped = set(prokka_skipped.keys()).union(
                set(skipped_faas.keys())
            )
            rows = []
            for bin_name in bins_skipped:
                row = []
                for c in results.columns:
                    if "score" in c or c == "entropy":
                        row.append(np.uint8(0.0))
                    else:
                        row.append("None")
                rows.append(row)
            skipped_df = pd.DataFrame(
                rows, index=bins_skipped, columns=results.columns
            )

            results = results.append(skipped_df)

    results_csv = output_dir / Path("results.csv")
    results.to_csv(results_csv, index_label="BIN/genome")

    print("**Deep learning predictions have finished. "
        "Results are in file {}.\n"
        "**Thank you for using vHULK".format(results_csv))

if __name__ == "__main__":
    main()
