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
from tensorflow.keras.models import load_model
from scipy.special import entr


def parse_arguments():

    parser = argparse.ArgumentParser(
        description="Predict phage draft genomes in metagenomic bins."
    )

    # Rearrange arguments to proper groups
    # to show up better in the help
    optionalArgs = parser._action_groups.pop()
    optionalArgs.title = "Optional arguments"

    # Required arguments
    requiredArgs = parser.add_argument_group("Required arguments")
    requiredArgs.add_argument(
        "-i",
        "--input_dir",
        required=True,
        type=lambda p: Path(p).resolve(strict=True),
        dest="input_dir",
        help="Path to a folder containing metagenomic bins in .fa or .fasta" "format",
    )

    requiredArgs.add_argument(
        "-o",
        "--output-dir",
        required=True,
        type=lambda p: Path(p).resolve(),
        dest="output_dir",
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
        default=Path("./models").resolve(),
        help="Path to directory where all models are stored.",
    )

    optionalArgs.add_argument(
        "-f",
        "--files-dir",
        required=False,
        dest="files_dir",
        default=Path("./files").resolve(),
        help="Files directory provided with vHULK",
    )

    parser._action_groups.append(optionalArgs)

    return parser.parse_args()


# Function declarations
# Get prefix from bins

# Run prokka
def run_prokka(fasta_in, output_dir, threads):
    # Check the fasta format
    fasta_suffix = fasta_in.suffix
    out_prefix = fasta_in.name.replace(fasta_suffix, "")
    genome_dir = output_dir / Path(out_prefix)
    # Filehandle where the output of prokka will be saved
    # output_prokka = open(str(prefix)+'prokka.output', mode='w')
    # Full command line for prokka
    command_line = (
        "prokka --kingdom Viruses --centre X --compliant "
        "--gcode 11 --cpus {} --force --quiet --prefix prokka_results_{} "
        "--fast --norrna --notrna --outdir {} "
        "--cdsrnaolap --noanno {}".format(threads, out_prefix, genome_dir, fasta_in)
    )

    return_code = subprocess.run(command_line, shell=True)

    return_code.check_returncode()


def run_hmmscan(fasta_in, output_dir, models_dir, threads):
    # Get the base name of the input file
    prefix = fasta_in.name.replace(fasta_in.suffix, "")
    basename = prefix.replace("prokka_results_", "")

    # Construct the filenames and paths
    ## for hmm stdout
    base_hmm_out = "_".join([str(basename), "hmmscan.out"])
    path_hmm_out = output_dir / Path(base_hmm_out)

    ## and hmm table output
    base_hmm_tblout = "_".join([str(basename), "hmmscan.tbl"])
    path_hmm_tblout = output_dir / Path(base_hmm_tblout)

    command_line = "hmmscan -o {} --cpu {} --tblout {} --noali {} " "{}".format(
        path_hmm_out, threads, path_hmm_tblout, models_dir, fasta_in
    )

    return_code = subprocess.run(command_line, shell=True)

    return_code.check_returncode()


def construct_gene_scores_matrix(hmmtable):
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


def predict_genus_relu(sliced_array, relu_genus_model, genus_hosts):
    pred_genus_relu = relu_genus_model.predict(sliced_array)
    idx_pred_genus_relu = np.argmax(pred_genus_relu)
    if not pred_genus_relu.any():
        name_genus_relu = "None"
        score_genus_relu = 0
    else:
        name_genus_relu = genus_hosts[idx_pred_genus_relu]
        score_genus_relu = pred_genus_relu[0][idx_pred_genus_relu]

    return name_genus_relu, score_genus_relu


def predict_genus_softmax(sliced_array, softmax_genus_model, genus_hosts):
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


def predict_species_relu(sliced_array, relu_species_model, species_hosts):
    pred_species_relu = relu_species_model.predict(sliced_array)
    idx_pred_species_relu = np.argmax(pred_species_relu)
    if not pred_species_relu.any():
        name_species_relu = "None"
        score_species_relu = 0
    else:
        name_species_relu = species_hosts[idx_pred_species_relu]
        score_species_relu = pred_species_relu[0][idx_pred_species_relu]

    return name_species_relu, score_species_relu


def predict_species_softmax(sliced_array, softmax_species_model, species_hosts):
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
    relu_genus_model,
    softmax_genus_model,
    relu_species_model,
    softmax_species_model,
    genus_hosts,
    species_hosts,
):

    predictions = {}

    for i in range(0, len(scores_array)):
        bin_name = names_array[i]
        bin_array = np.array([scores_array[i]])
        name_genus_relu, score_genus_relu = predict_genus_relu(
            bin_array, relu_genus_model, genus_hosts
        )

        name_genus_sm, score_genus_sm, entropy_genus_sm = predict_genus_softmax(
            bin_array, softmax_genus_model, genus_hosts
        )

        name_species_relu, score_species_relu = predict_species_relu(
            bin_array, relu_species_model, species_hosts
        )

        name_species_sm, score_species_sm = predict_species_softmax(
            bin_array, softmax_species_model, species_hosts
        )

        # Decision Tree
        final_decision = "None"
        # Relu sp
        if score_species_relu > 0.9:
            final_decision = name_species_relu
        # SM sp
        if (score_species_sm > 0.6) and (name_species_sm != final_decision):
            final_decision = name_species_sm
        # Coudn't predict species
        if final_decision == "None":
            # Put here sm sp
            if score_species_sm > 0.6:
                final_decision = name_species_sm
                # relu genus
                if score_genus_relu >= 0.7:
                    final_decision = name_genus_relu
                # sm genus
                if (score_genus_sm >= 0.5) and (name_genus_sm != final_decision):
                    final_decision = name_genus_sm
            else:
                # relu genus
                if score_genus_relu >= 0.9:
                    final_decision = name_genus_relu
                # sm genus
                if (score_genus_sm >= 0.4) and (name_genus_sm != final_decision):
                    final_decision = name_genus_sm
        # Predicted species.
        # Verify if genus is the same
        else:
            if re.search(name_genus_relu, final_decision) or re.search(
                name_genus_sm, final_decision
            ):
                pass
            else:
                # relu genus
                if score_genus_relu >= 0.9:
                    final_decision = name_genus_relu
                # sm genus
                if (score_genus_sm >= 0.5) and (name_genus_sm != final_decision):
                    final_decision = name_genus_sm

        predictions[bin_name] = {
            "pred_genus_relu": name_genus_relu,
            "score_genus_relu": score_genus_relu,
            "pred_genus_softmax": name_genus_sm,
            "score_genus_softmax": score_genus_sm,
            "pred_species_relu": name_species_relu,
            "score_species_relu": score_species_relu,
            "pred_species_softmax": name_species_sm,
            "score_species_softmax": score_species_sm,
            "final_prediction": final_decision,
            "entropy": entropy_genus_sm,
        }

    return predictions


####
#### Main code
####


def main():
    args = parse_arguments()
    # Greeting message
    print("\n**Welcome v.HULK, a toolkit for phage host prediction!\n")
    now = datetime.datetime.now()
    print(now.strftime("**%Y-%m-%d %H:%M:%S"))
    # Verify databases
    vog_profiles = args.models_dir / Path("all_vogs_hmm_profiles_feb2018.hmm")

    if not vog_profiles.exists():
        print(
            "**Your database and models are not set. "
            "Please, run: python download_and_set_models.py \n"
        )
        sys.exit(1)

    ## Important variables

    input_dir = args.input_dir
    output_dir = args.output_dir
    models_dir = args.models_dir
    files_dir = args.files_dir
    threads = args.threads

    # Get a list of all input fastas from the dir
    list_bins = []
    for entry in input_dir.iterdir():
        if entry.is_file() and entry.suffix.lower() in [".fa", ".fasta"]:
            list_bins.append(entry)

    # Check if there is some input
    if len(list_bins) == 0:
        print(
            "**Input folder appears to be empty or does not"
            "contain expected bin files (suffixed with 'fa' or '.fasta')."
            "Exiting...\n"
        )
        sys.exit(1)
    else:

        print(
            "**Arguments are OK. Checked the input folder "
            "and found {} genome files.\n".format(len(list_bins))
        )

    # Check output dir and create it if needed
    if not output_dir.is_dir():
        output_dir.mkdir(parents=True, exist_ok=True)

    ######
    ## PROKKA
    ######
    count_prokka = 0
    prokka_skipped = []
    prokka_dir = output_dir / Path("prokka")

    print("**Prokka has started, this may take a while. Be patient.\n")
    for bin_fasta in list_bins:
        len_bin = 0
        for record in SeqIO.parse(bin_fasta, "fasta"):
            len_bin += len(record.seq)
        if len_bin < 5000:
            print(
                "**v.HULK has found a genome or bin, which is too short to "
                "code proteins (< 5000 bp). As CDSs are an import feature for "
                "v.HULK, we will be skipping this: " + bin_fasta.name
            )
            prokka_skipped.append(bin_fasta.name)
            continue

        run_prokka(bin_fasta, prokka_dir, threads)

        count_prokka += 1
        if count_prokka % 10 == 0:
            print("**Done with {} genomes...".format(count_prokka))

    print("**Successfully run: {}".format(count_prokka))

    if len(prokka_skipped) != 0:
        print("**Skipped : {}".format(len(prokka_skipped)))

    now = datetime.datetime.now()
    print(now.strftime("**%Y-%m-%d %H:%M:%S"))
    #####
    ## HMM SEARCHES
    #####

    print("**Starting HMM scan, this may take awhile. Be patient.\n")

    hmmscan_dir = output_dir / Path("hmmscan")
    if not hmmscan_dir.is_dir():
        hmmscan_dir.mkdir(parents=True, exist_ok=True)

    # Filter out empty protein faas
    # This might occur when prokka doesn't call any genes
    # TO DO
    # Check with prokka if this is desired behavior
    valid_faas, skipped_faas = {}, {}
    for faa in prokka_dir.glob("**/*.faa"):
        suffix = faa.suffix
        bin_name = faa.name.replace(suffix, "").replace("prokka_results_", "")
        if faa.stat().st_size != 0:
            valid_faas[bin_name] = faa
        else:
            skipped_faas[bin_name] = faa

    count_hmms = 0
    for faa in valid_faas.values():
        run_hmmscan(faa, hmmscan_dir, vog_profiles, threads)
        count_hmms += 1
        print("**Done with {} / {} HMMs\r".format(count_hmms, len(valid_faas)), end="")
    else:
        print("\n**Done with HMMscan!")

    now = datetime.datetime.now()
    print(now.strftime("**%Y-%m-%d %H:%M:%S"))

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
            index=gene_scores.keys(), columns=dic_vogs_headers.keys(), dtype=float
        )

        dic_matrices_by_genome[entry].fillna(value=np.float32(0.0), inplace=True)
        for gene in gene_scores:
            for each_match in gene_scores[gene]:
                if each_match[1] > 1:
                    each_match[1] = 1
                dic_matrices_by_genome[entry][each_match[0]][gene] = np.float32(
                    1.0 - np.float32(each_match[1])
                )

    list_condensed_matrices = []
    bin_names = []
    for matrix in dic_matrices_by_genome:
        temp = list(dic_matrices_by_genome[matrix].sum(axis=0, skipna=True))
        bin_names.append(matrix)
        list_condensed_matrices.append(temp)

    mat_array = np.array(list_condensed_matrices)

    #################
    ## Predictions ##
    #################
    print("**Starting deep learning predictions")
    print("**Loading models...")

    # genus relu
    genus_relu_h5 = models_dir / Path("model_genus_total_fixed_relu_08mar_2020.h5")

    model_genus_relu = load_model(
        str(genus_relu_h5), custom_objects={"LeakyReLU": LeakyReLU, "ReLU": ReLU}
    )
    # genus softmax
    genus_sm_h5 = models_dir / Path("model_genus_total_fixed_softmax_01mar_2020.h5")

    model_genus_sm = load_model(
        str(genus_sm_h5), custom_objects={"LeakyReLU": LeakyReLU, "ReLU": ReLU}
    )
    # species relu
    species_relu_h5 = models_dir / Path("model_species_total_fixed_relu_08mar_2020.h5")
    model_species_relu = load_model(
        str(species_relu_h5), custom_objects={"LeakyReLU": LeakyReLU, "ReLU": ReLU}
    )
    # species softmax
    species_sm_h5 = models_dir / Path("model_species_total_fixed_softmax_01mar_2020.h5")
    model_species_sm = load_model(
        str(species_sm_h5), custom_objects={"LeakyReLU": LeakyReLU, "ReLU": ReLU}
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
        model_genus_relu,
        model_genus_sm,
        model_species_relu,
        model_species_sm,
        genus_hosts,
        species_hosts,
    )

    # Put it in a dataframe and write it
    results_df = pd.DataFrame.from_dict(a, orient="index")
    results_csv = output_dir / Path("results.csv")
    results_df.to_csv(results_csv, index_label="BIN/genome")

    print(
        "\n**Deep learning predictions have finished. "
        "Results are in file {}.\n"
        "**Thank you for using v.HULK".format(results_csv)
    )


if __name__ == "__main__":
    main()
