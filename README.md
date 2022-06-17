<p align="center"> <img src="logo_vHULK.png" height="200" alt="vHULK" /> </p>

# vHULK

**V**iral **H**ost **U**nvei**L**ing **K**it - A toolkit for phage host 
prediction

### Pre-print publication

A preview of vHULK's publication is available in 
[BioRxiv](https://doi.org/10.1101/2020.12.06.413476) 
as the original manuscript is under review.

## Release Note June, 17th 2022

vHULK is going through updates. Consider only the scripts **download_and_set_models.py** and **training.py**.

Shortly we will reimplement the tool vHULK.

## Release Note March, 27th 2021

We implemented some changes to the tool following suggestions provided by Reviewers. These include:
* L2 regularization to weights of the Artificial Neural Network
* Models based on ReLU activation functions were removed
* Models now are included in this git repository and will be downloaded along with source code (cloning)
* The decision tree used to make an unified final prediction was simplified
* Other minor changes

## Phage Host Prediction using high level features and neural networks

Metagenomics and sequencing techniques have greatly improved in these last 
five years and, as a consequence, the amount of data from microbial communities
is astronomic. An import part of the microbial community are phages, which have
their own ecological roles in the environment. Besides that, they have also 
been given a possible human relevant (clinical) role as terminators of multidrug
resistant bacterial infections. A lot of basic research still need to be done 
in the **Phage therapy** field, and part of this research involves gathering 
knowledge from new phages present in the environment as well as about their 
relationship with clinical relevant bacterial pathogens.

Having this scenario in mind, we have developed vHULK. A user-friendly tool 
for prediction of phage hosts given their complete or partial genome in FASTA 
format. Our tool outputs an ensemble prediction at the genus or species level 
based on scores of four different neural network models. Each model was trained
with more than 4,000 genomes whose phage-host relationship was known. v.HULK 
also outputs a mesure of entropy for each final prediction, which we have 
demonstrated to be correlated with prediction's accuracy. The user might 
understand this value as additional information of how certain v.HULK is about
a particular prediction. We also suspect that phages with higher entropy values
may have a broad host-range. But that hypothesis is to be tested later. 
**Accuracy results in test datasets were >99% for predictions at the genus 
level and >98% at the species level**. vHULK currently supports predictions for
52 different prokaryotic host species and 61 different genera.

Target species
```
Acinetobacter_baumannii, Aeromonas_salmonicida, Arthrobacter_globiformis, Bacillus_cereus, Bacillus_megaterium, Bacillus_pumilus, Bacillus_subtilis, Bacillus_thuringiensis, Brucella_abortus, Campylobacter_jejuni, Caulobacter_crescentus, Cellulophaga_baltica, Citrobacter_freundii, Clostridium_difficile, Clostridium_perfringens, Cronobacter_sakazakii, Enterococcus_faecalis, Erwinia_amylovora, Escherichia_coli, Gordonia_terrae, Helicobacter_pylori, Klebsiella_pneumoniae, Lactococcus_lactis, Listeria_monocytogenes, Microbacterium_foliorum, Moraxella_catarrhalis, Mycobacterium_smegmatis, Paenibacillus_larvae, Pectobacterium_atrosepticum, Propionibacterium_acnes, Proteus_mirabilis, Pseudomonas_aeruginosa, Pseudomonas_fluorescens, Pseudomonas_syringae, Ralstonia_solanacearum, Rhodococcus_erythropolis, Salmonella_enterica, Salmonella_typhimurium, Shigella_flexneri, Shigella_sonnei, Staphylococcus_aureus, Staphylococcus_epidermidis, Streptococcus_pneumoniae, Streptococcus_thermophilus, Streptomyces_griseus, Streptomyces_lividans, Streptomyces_venezuelae, Sulfolobus_islandicus, Vibrio_anguillarum, Vibrio_cholerae, Vibrio_parahaemolyticus, Yersinia_enterocolitica
```
Target genera
```
Achromobacter, Acidianus, Acinetobacter, Aeromonas, Arthrobacter, Bacillus, Brevibacillus, Brucella, Burkholderia, Campylobacter, Caulobacter, Cellulophaga, Citrobacter, Clostridium, Corynebacterium, Cronobacter, Dickeya, Enterobacter, Enterococcus, Erwinia, Escherichia, Flavobacterium, Gordonia, Haloarcula, Halorubrum, Helicobacter, Klebsiella, Lactobacillus, Lactococcus, Leuconostoc, Listeria, Mannheimia, Microbacterium, Moraxella, Mycobacterium, Paenibacillus, Pectobacterium, Prochlorococcus, Propionibacterium, Proteus, Pseudoalteromonas, Pseudomonas, Ralstonia, Rhizobium, Rhodococcus, Ruegeria, Salmonella, Serratia, Shigella, Sinorhizobium, Staphylococcus, Stenotrophomonas, Streptococcus, Streptomyces, Sulfolobus, Synechococcus, Thermus, Vibrio, Xanthomonas, Xylella, Yersinia
```

### Scripts

Main script:
   * **vHULK.py** - Deep learning host prediction from phage genomes
  
Auxiliary script:
   * **download_and_set_models.py** - Used only at the first use for 
   downloading and setting models
   
### Dependencies

All scripts from this project were coded in 
[Python 3.7](https://www.python.org/). So, first of all, make sure you have it
installed and updated.

vHULK's main script (`vHULK.py`) requires `Prokka` and `HMMER` tools as 
dependencies. 

* [Prokka](https://github.com/tseemann/prokka) - Rapid Prokaryotic genome annotation.
* [HMMER Tools](http://www.hmmer.org/) - Biosequence analysis using profile hidden Markov models

These python libraries are required:

* [Numpy](http://www.numpy.org/), [Pandas](https://pandas.pydata.org/), 
* [Scipy](https://www.scipy.org/) - Efficiently handling arrays and scientific computing for Python
* [Biopython](http://biopython.org/) - Handling biological sequences
* [Tensorflow](https://www.tensorflow.org/) - Google's Deep Neural Networks libraries 

To install these dependencies, there are usually two ways: 
[pip](https://pypi.org/project/pip/) or 
[conda](https://www.anaconda.com/products/individual).  
We strongly recommend creating a specific conda environment containing the installed libraries and tools.

> The `$` signifies the command-line prompt

Install all requirements using these commands:
```
$ conda create -n vHULK -c conda-forge -c bioconda -c defaults prokka
$ conda install -n vHULK -c bioconda hmmer
$ conda install -n vHULK -c bioconda -c anaconda numpy pandas scipy biopython tensorflow=2.4.1
$ conda activate vHULK 
```

> Note
>
> Some people have been facing issues with newer versions of Prokka or TensorFlow. 
> To avoid such problems, we froze a working version with TensorFlow 2.4.1 we know
> for sure it works with vHULK.

If all went well your command-line prompt will now be prefixed with the 
environment name you gave above.

```
(vHULK)$ python -V
Python 3.7.X
```

## Installation

Getting vHULK ready to run is as simple as cloning this Github project or 
download and extract it to a directory inside your computer:

```
$ git clone https://github.com/LaboratorioBioinformatica/vHULK
```

### Download models

The helper script `download_and_set_models.py` can be used to download all 
data dependencies required for vHULK to run.

This is required only once and it is simple. Just run:
```
(vHULK)$ python download_and_set_models.py
```
This will download HMM models to a directory called `models`, which stores all necessary 
data. You can move the whole directory in a location of your preference and
point to it with the `-m` option of `vHULK.py` when running vHULK. 
By default, vHULK searches for it in the same path where the executable is 
located.


## Usage

The `vHULK.py` is an executable script. Its location can be included in your 
`$PATH` environmet variable if you so desire.

For example, assuming you have cloned this repo in `/home/user/tools/vHULK`,
you can prepend this location to `$PATH`
```
$ export PATH="/home/user/tools/vHULK:$PATH"
```

This allows you to invoke `vHULK.py` from any location on your system.

Alternatively, all invocations must point to its full location e.g. 
```
(vHULK)$ /full/path/to/vHULK.py -h
```

Here, it is assumed tha it is available in your `$PATH`.

To list all options available run
```
(vHULK)$ vHULK.py -h
usage: vHULK.py [-h] -i INPUT_DIR -o OUTPUT_DIR [-t THREADS] [-m MODELS_DIR]
                [-f FILES_DIR] [--all] [-v]

Predict phage draft genomes in metagenomic bins.

Required arguments:
  -i INPUT_DIR, --input-dir INPUT_DIR
                        Path to a folder containing metagenomic bins in .fa or
                        .fasta format (default: )
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Location to store results in. It is created if it
                        doesn't exist (default: )

Optional arguments:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        Number of CPU threads to be used by Prokka and hmmscan
                        (default: 1)
  -m MODELS_DIR, --models-dir MODELS_DIR
                        Path to directory where all models are stored.
                        (default: /home/user/tools/vHULK/models)
  -f FILES_DIR, --files-dir FILES_DIR
                        Files directory provided with vHULK (default:
                        /home/user/tools/Projects/vHULK/files)
  --all                 Write predictions for all input bins/genomes, even if
                        they were skipped (size filtered or hmmscan failed)
                        (default: False)
  -v, --version         show program's version number and exit
```

### Running the example datasets

We provide a folder with example datasets containing mocking bins of RefSeq 
viral and bacterial genomes.  

To try these examples, run:

```
(vHULK)$ vHULK.py -i test_input -o test_output -t 4
```

It should take about 2 minutes to generate your nice and accurate predictions.

## Input

vHULK is ready to accept whole or partial phage genomes. Keep in mind that vHULK's 
predictions are based in high level annotated features, i.e., features that 
depend on gene annotation. So, very small contigs with few or no entire genes 
will not provide features for vHULK to work with. In general, if you have a 
partial genome larger than 10 kbp you should be fine. Currently, vHULK skips 
all genomes or bins that are smaller than 5000bp at runtime.

Only FASTA nucleotide files are accepted, and vHULK will understand each 
individual file as one phage genome. Therefore, single-sequence FASTA files 
will be understood as whole or partial phage genomes and multiFASTA as a 
metagenomic BIN. Our tests indicate that prediction's accuracy is not affected 
in fragmented genomes (BINs).

## Output

vHULK's main output file is `results.csv`, which will be inside the output 
directory you provided with the `-o` option, after execution is finished. 
If more than one fasta file is inside the input_folder, vHULK will understand 
that there are more than one bin or genome for prediction. Therefore, each line
of the CSV file will correspond to a bin or genome. 
Identifiers will be in the first column.

Output example:
```
BIN/genome,pred_genus_softmax,score_genus_softmax,pred_species_softmax,score_species_softmax,final_prediction,entropy
MK524521_mycobacterium_smegmatis,Mycobacterium,0.9322084188461304,Mycobacterium_smegmatis,0.9966434240341187,Mycobacterium_smegmatis,0.2768826484680176
MK570225_enterococcus_faecalis,Enterococcus,0.9962823987007141,Enterococcus_faecalis,0.9984235763549805,Enterococcus_faecalis,0.033856652677059174
MK801680_BIN_staphylococcus_aureus,Staphylococcus,0.9952540397644043,Staphylococcus_aureus,0.9958204030990601,Staphylococcus_aureus,0.040766458958387375
MN41915_vibrio_cholerae,Vibrio,0.9746623635292053,Vibrio_cholerae,0.9812621474266052,Vibrio_cholerae,0.17507237195968628
MN689520_lactococcus_lactis,Lactococcus,0.981468915939331,Lactococcus_lactis,0.9516076445579529,Lactococcus_lactis,0.14059613645076752
ZC4_composting_unknown,Streptomyces,0.8378896713256836,Streptomyces_lividans,0.6341040134429932,Streptomyces_lividans,0.955293595790863
```

Note that there is a header in the output file, and it is self explanatory. 

Nonetheless:  

* 1st column: ID of the submitted sequence  
* 2nd to 5th columns: Host prediction and respective score for each one of the
two vHULK's internal models  
* 6th column: Unified host prediction generated by an heuristic decision tree
based on score values  
* 7th column: Entropy value for predictions of the Softmax model at the genus 
level. It can be used as a proxy of confidence for a particular prediction, 
i.e., values close to 0 correlate with higher accuracies and values close to 4 
otherwise.

The output directory also contains two more subdirectories holding the 
intermediate results of prokka and hmmscan. These are named accordingly.

