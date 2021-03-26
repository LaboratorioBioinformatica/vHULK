#!/usr/bin/env python
# coding: utf-8
# Edited May, 27th 2020
## This is vHULK: viral Host Unveiling Kit
# Developed by Deyvid Amgarten and Bruno Iha
# Creative commons

# Import required Python modules
import numpy as np
import pandas as pd
from Bio import SeqIO, SearchIO
import re
import sys
import os

# Set logging level fro TensorFlow
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"

from pathlib import Path
import subprocess
import datetime
import argparse

# Suppress warnings from Tensorflow
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.simplefilter(action="ignore", category=FutureWarning)

# Deep Learning
from tensorflow.keras.layers import Dense, Activation, LeakyReLU, ReLU
from scipy.special import entr

### Adicionadas por Bruno em 03/08/2021 para fazer Load dos Modelos
from tensorflow.keras.models import Sequential, load_model
from tensorflow.keras import optimizers
from tensorflow.keras.regularizers import l2
import tensorflow.keras.utils as utils
from tensorflow.keras.utils import Sequence
from tensorflow.keras import optimizers, initializers
from tensorflow.keras.optimizers import Adam
###

__version__ = "1.0.0"

## REQUIRED FILES
# models
MODEL_FILES = [
    "GenusModel03Mar2021",
    "SpeciesModel03Mar2021"
]

# vogs hmms
VOG_PROFILES = "all_vogs_hmm_profiles_feb2018.hmm"

# host annotation
HOST_FILES = [
    "list_hosts_genus.txt",
    "list_hosts_species.txt",
    "VOGs_header.txt",
]

use_downloads_msg = (
    "Please run \n`python download_and_set_models.py` \n"
    "to make sure you have all input data required for v.HULK\n"
)

# Full path to this script
here = Path(__file__).resolve().parent


def parse_arguments():

    parser = argparse.ArgumentParser(
        prog="vHULK.py",
        description="Predict phage draft genomes in metagenomic bins.",
        # Display default values when printing help
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Rearrange arguments to proper groups
    # to show up better in the help
    optionalArgs = parser._action_groups.pop()
    optionalArgs.title = "Optional arguments"

    # Required arguments
    requiredArgs = parser.add_argument_group("Required arguments")
    requiredArgs.add_argument(
        "-i",
        "--input-dir",
        required=True,
        type=lambda p: Path(p).resolve(strict=True),
        dest="input_dir",
        default="",
        help="Path to a folder containing metagenomic bins in .fa or .fasta "
        "format",
    )

    requiredArgs.add_argument(
        "-o",
        "--output-dir",
        required=True,
        type=lambda p: Path(p).resolve(),
        dest="output_dir",
        default="",
        help="Location to store results in. It is created if it doesn't exist",
    )

    # Optional args
    optionalArgs.add_argument(
        "-t",
        "--threads",
        dest="threads",
        required=False,
        default=1,
        type=int,
        help="Number of CPU threads to be used by Prokka and hmmscan",
    )

    optionalArgs.add_argument(
        "-m",
        "--models-dir",
        required=False,
        dest="models_dir",
        type=lambda p: Path(p).resolve(),
        default=here / Path("models"),
        help="Path to directory where all models are stored.",
    )

    optionalArgs.add_argument(
        "-f",
        "--files-dir",
        required=False,
        dest="files_dir",
        default=here / Path("files"),
        help="Files directory provided with vHULK",
    )

    optionalArgs.add_argument(
        "--all",
        required=False,
        dest="write_all",
        action="store_true",
        help="Write predictions for all input bins/genomes, even if they "
        "were skipped (size filtered or hmmscan failed)",
    )

    optionalArgs.add_argument(
        "-v",
        "--version",
        action="version",
        version=__version__,
    )

    parser._action_groups.append(optionalArgs)

    return parser.parse_args()


def get_bin_name(fasta_p):
    """
    Strip relevant pre and suffixes from a given fasta name

    Arguments:
        fasta_p: pathlib.Path instance: The path to the fasta file
    Return:
        bin_name: str: The name of the bin as string
    """
    bin_name = fasta_p.name
    fa_suffix = fasta_p.suffix
    prokka_prefix = "prokka_results_"
    # Check if it is prefixed with prokka_results
    if bin_name.startswith(prokka_prefix):
        bin_name = fasta_p.name.replace(prokka_prefix, "")
    # Remove the .fa[a|sta] suffix
    bin_name = bin_name.replace(fa_suffix, "")
    return bin_name


# TO DO
# Use logging
def print_now():
    """
    Prints current date and time
    """
    now = datetime.datetime.now()
    print(now.strftime("\n== %Y-%m-%d %H:%M:%S ==\n"))


def run_prokka(fasta_in, output_dir, threads):
    """
    Run prokka for a given fasta file.

    Raises CalledProcessError if command doesn't end succesfully

    Arguments:
        fasta_in: pathlib.Path instance: Path to fasta file
        output_dir: pathlib.Path instance: Path to directory where results will
            be stored
        threads: int: Number of cpus/threads for prokka to use
    Return:
        -
    """
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


def run_hmmscan(fasta_in, output_dir, vogs_hmms, threads):
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


def construct_gene_scores_matrix(hmmtable):
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
    dic_genes_scores = {}
    for gene in SearchIO.parse(hmmtable, "hmmer3-tab"):
        dic_genes_scores[gene.id] = []
        for hit in gene.hits:
            hit_info = [
                hit.id,
                np.float32(hit.evalue),
                np.float32(hit.bitscore),
                np.float32(hit.bias),
            ]
            dic_genes_scores[gene.id].append(hit_info)
    return dic_genes_scores

def predict_genus_softmax(sliced_array, softmax_genus_model, genus_hosts):
    """
    Predict genus based on the genus softmax model
    """
    pred_genus_sm = softmax_genus_model.predict(sliced_array)
    # Calculate entropy
    entropy_genus_sm = entr(pred_genus_sm).sum(axis=1)

    idx_pred_genus_sm = np.argmax(pred_genus_sm)

    if not pred_genus_sm.any():
        name_genus_sm = "None"
        score_genus_sm = 0
    else:
        name_genus_sm = genus_hosts[idx_pred_genus_sm]
        score_genus_sm = pred_genus_sm[0][idx_pred_genus_sm]

    return name_genus_sm, score_genus_sm, entropy_genus_sm[0]

def predict_species_softmax(
    sliced_array, softmax_species_model, species_hosts
):
    """
    Predict species based on the species softmax model
    """
    pred_species_sm = softmax_species_model.predict(sliced_array)
    idx_pred_species_sm = np.argmax(pred_species_sm)
    if not pred_species_sm.any():
        name_species_sm = "None"
        score_species_sm = 0
    else:
        name_species_sm = species_hosts[idx_pred_species_sm]
        score_species_sm = pred_species_sm[0][idx_pred_species_sm]

    return name_species_sm, score_species_sm

def predict(
    scores_array,
    names_array,
    softmax_genus_model,
    softmax_species_model,
    genus_hosts,
    species_hosts,
):
    """
    Predict genus and species for all input for all models
    """

    predictions = {}

    for i in range(0, len(scores_array)):
        bin_name = names_array[i]
        bin_array = np.array([scores_array[i]])
        (
            name_genus_sm,
            score_genus_sm,
            entropy_genus_sm,
        ) = predict_genus_softmax(bin_array, softmax_genus_model, genus_hosts)

        name_species_sm, score_species_sm = predict_species_softmax(
            bin_array, softmax_species_model, species_hosts
        )

        # Decision Tree
        final_decision = "None"
        # SM sp
        if (score_species_sm > 0.5) and (name_species_sm != final_decision):
            final_decision = name_species_sm
        # Coudn't predict species
        if final_decision == "None":
            # Put here sm sp
            if score_species_sm > 0.5:
                final_decision = name_species_sm
                # sm genus
                if (score_genus_sm >= 0.5) and (
                    name_genus_sm != final_decision
                ):
                    final_decision = name_genus_sm
            else:
                # sm genus
                if (score_genus_sm >= 0.5) and (
                    name_genus_sm != final_decision
                ):
                    final_decision = name_genus_sm
        # Predicted species.
        # Verify if genus is the same
        else:
            if re.search(name_genus_sm, final_decision):
                pass
            else:
                # sm genus
                if (score_genus_sm >= 0.5) and (
                    name_genus_sm != final_decision
                ):
                    final_decision = name_genus_sm

        predictions[bin_name] = {
            "pred_genus_softmax": name_genus_sm,
            "score_genus_softmax": score_genus_sm,
            "pred_species_softmax": name_species_sm,
            "score_species_softmax": score_species_sm,
            "final_prediction": final_decision,
            "entropy": entropy_genus_sm,
        }

    return predictions


##########
## MAIN ##
##########


def main():
    args = parse_arguments()

    # Greeting message
    print("\n**Welcome to vHULK, a toolkit for phage host prediction!")

    print_now()
    ## Important variables

    input_dir = args.input_dir
    output_dir = args.output_dir
    models_dir = args.models_dir
    files_dir = args.files_dir
    threads = args.threads

    if not models_dir.is_dir():
        print("ERROR :Missing models directory {}\n".format(models_dir))
        print(use_downloads_msg)
        sys.exit(1)

    # Verify databases
    files_ok = []
    for m in MODEL_FILES:
        model_path = models_dir / Path(m)
        is_there = model_path.exists()
        files_ok.append(is_there)
        if not is_there:
            print("ERROR: Missing {} from {}".format(model_path, models_dir))

    vog_profiles = models_dir / Path(VOG_PROFILES)
    vogs_there = vog_profiles.exists()
    files_ok.append(vogs_there)
    if not vogs_there:
        print("ERROR: Missing {} from {}".format(vog_profiles, models_dir))

    for h in HOST_FILES:
        txt_fp = files_dir / Path(h)
        h_there = txt_fp.exists()
        files_ok.append(h_there)
        if not h_there:
            print("ERROR: Missing {} from {}".format(txt_fp, files_dir))

    if not all(files_ok):
        print(use_downloads_msg)
        sys.exit(1)

    # Get a list of all input fastas from the dir
    list_bins = []
    for entry in input_dir.iterdir():
        if entry.is_file() and entry.suffix.lower() in [".fa", ".fasta"]:
            list_bins.append(entry)

    # Check if there is some input
    if len(list_bins) == 0:
        print(
            "**ERROR: Input folder appears to be empty or does not"
            "contain expected bin files (suffixed with 'fa' or '.fasta')."
        )
        sys.exit(1)
    else:

        print(
            "**Arguments are OK. Checked the input folder "
            "and found {} genome files.".format(len(list_bins))
        )

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
            print(
                "**vHULK has found a genome or bin, which is too short to "
                "code proteins (< 5000 bp). As CDSs are an import feature for "
                "vHULK, we will be skipping this: " + bin_fasta.name
            )
            prokka_skipped[bin_name] = bin_fasta
            continue

        run_prokka(bin_fasta, prokka_dir, threads)

        count_prokka += 1
        if count_prokka % 10 == 0:
            print("**Done with {} genomes...".format(count_prokka))

    print("\n**PROKKA finished with no errors")
    print(
        "**{:>{width}} : {}".format("Successfully run", count_prokka, width=20)
    )

    if len(prokka_skipped) != 0:
        print(
            "**{:>{width}} : {}".format(
                "Skipped", len(prokka_skipped), width=20
            )
        )

    print_now()

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
        print(
            "**ERROR: No valid protein fastas were found in {}".format(
                prokka_dir
            )
        )
        sys.exit(1)

    ##################
    ## HMM SEARCHES ##
    ##################

    print("**Starting HMM scan, this may take a while. Be patient.")

    hmmscan_dir = output_dir / Path("hmmscan")
    if not hmmscan_dir.is_dir():
        hmmscan_dir.mkdir(parents=True, exist_ok=True)

    count_hmms = 0
    for faa in valid_faas.values():
        run_hmmscan(faa, hmmscan_dir, vog_profiles, threads)
        count_hmms += 1
        print(
            "**Done with {} / {} HMMs\r".format(count_hmms, len(valid_faas)),
            end="",
        )
    else:
        print("\n**Done with HMMscan!")

    print_now()

    # Initialize a dictionary that holds some info
    dic_vogs_headers = {}
    headers_fp = files_dir / Path("VOGs_header.txt")
    with open(headers_fp, "r") as fin:
        for line in fin:
            vog = line.rstrip()
            dic_vogs_headers[vog] = np.float32(0.0)

    print("**Parsing HMM results")
    # PARSE HMM RESULTS
    dic_matrices_by_genome = {}
    for entry in valid_faas:
        hmmtbl_str = "_".join([entry, "hmmscan.tbl"])
        hmmtable = hmmscan_dir / Path(hmmtbl_str)
        gene_scores = construct_gene_scores_matrix(hmmtable)
        dic_matrices_by_genome[entry] = pd.DataFrame(
            index=gene_scores.keys(),
            columns=dic_vogs_headers.keys(),
            dtype=float,
        )

        dic_matrices_by_genome[entry].fillna(
            value=np.float32(0.0), inplace=True
        )
        for gene in gene_scores:
            for each_match in gene_scores[gene]:
                if each_match[1] > 1:
                    each_match[1] = 1
                dic_matrices_by_genome[entry][each_match[0]][
                    gene
                ] = np.float32(1.0 - np.float32(each_match[1]))

    list_condensed_matrices = []
    bin_names = []
    for matrix in dic_matrices_by_genome:
        temp = list(dic_matrices_by_genome[matrix].sum(axis=0, skipna=True))
        bin_names.append(matrix)
        list_condensed_matrices.append(temp)

    mat_array = np.array(list_condensed_matrices)
    print_now()
    
    #################
    ## PREDICTIONS ##
    #################
    print("**Starting deep learning predictions")

    # genus softmax
    genus_sm_h5 = models_dir / Path(
        "GenusModel03Mar2021"
    )

    model_genus_sm = load_model(
        str(genus_sm_h5)
    )

    # species softmax
    species_sm_h5 = models_dir / Path(
        "SpeciesModel03Mar2021"
    )
    model_species_sm = load_model(
        str(species_sm_h5)
    )
    # Read in the hosts lists
    ## Genus hosts
    genus_hosts_fp = files_dir / Path("list_hosts_genus.txt")
    with open(genus_hosts_fp, "r") as fin:
        genus_hosts = [line.rstrip() for line in fin]

    ## Species hosts
    species_hosts_fp = files_dir / Path("list_hosts_species.txt")
    with open(species_hosts_fp, "r") as fin:
        species_hosts = [line.rstrip() for line in fin]

    # Predict them all and store in a dictionary
    a = predict(
        mat_array,
        bin_names,
        model_genus_sm,
        model_species_sm,
        genus_hosts,
        species_hosts,
    )

    # Put it in a dataframe and write it
    results_df = pd.DataFrame.from_dict(a, orient="index")

    if args.write_all:
        if len(prokka_skipped) != 0 or len(skipped_faas) != 0:
            bins_skipped = set(prokka_skipped.keys()).union(
                set(skipped_faas.keys())
            )
            rows = []
            for bin_name in bins_skipped:
                row = []
                for c in results_df.columns:
                    if "score" in c or c == "entropy":
                        row.append(np.float32(0.0))
                    else:
                        row.append("None")
                rows.append(row)
            skipped_df = pd.DataFrame(
                rows, index=bins_skipped, columns=results_df.columns
            )

            results_df = results_df.append(skipped_df)

    results_csv = output_dir / Path("results.csv")
    results_df.to_csv(results_csv, index_label="BIN/genome")
    print_now()
    print(
        "**Deep learning predictions have finished. "
        "Results are in file {}.\n"
        "**Thank you for using vHULK".format(results_csv)
    )


if __name__ == "__main__":
    main()