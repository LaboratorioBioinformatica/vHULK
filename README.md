<p align="center"> <img src="logo_vHULK.png" height="110" alt="vHULK" /> </p>
# v.HULK
Viral Host UnveiLing Kit - A tool-kit for phage host prediction

### Phage Host Prediction using high level features and neural networks
Metagenomics and sequencing techiniques have greatly improved in these last five years and, as a consequence, the amount of data from microbial communities is astronomic. An import part of the microbial community are phages, which have their own ecological roles in the environment. Besides that, they have been also given a possible human relevant (clinical) role as terminators of multidrug resistant bacterial infections. A lot of basic research still need to be done in the **Phage therapy** field, and part of this research involves gathering knowledge from new phages present in the environment as well as about their relationship with clinical relevant bacterial pathogens.

We have developed v.HULK having this scenario in mind. An user-friendly tool for prediction of phage hosts given their complete or partial genome in FASTA format. Our tool outputs an ensamble prediction at the genus or species level based on scores of four different neural network models. Each model was trained with more than 4,000 genomes whose phage-host relatonship was known. v.HULK outputs a mesure of entropy, which we have demonstrated to be correlated with prediction's accuracy. The user might understand this value as additional information of how certain v.HULK is about a particular prediction. We also suspect that phage with higher entropy values may have a broad host-range. But that hypothesis is to be tested later. **Accuracy results in test datasets were >99% for predictions at the genus level and >98% at the species level**. vHULK currently supports predictions for 52 different prokaryotic host species and 62 different genera.

Target species
```
Acinetobacter_baumannii, Aeromonas_salmonicida, Arthrobacter_globiformis, Bacillus_cereus, Bacillus_megaterium, Bacillus_pumilus, Bacillus_subtilis, Bacillus_thuringiensis, Brucella_abortus, Campylobacter_jejuni, Caulobacter_crescentus, Cellulophaga_baltica, Citrobacter_freundii, Clostridium_difficile, Clostridium_perfringens, Cronobacter_sakazakii, Enterococcus_faecalis, Erwinia_amylovora, Escherichia_coli, Gordonia_terrae, Helicobacter_pylori, Klebsiella_pneumoniae, Lactococcus_lactis, Listeria_monocytogenes, Microbacterium_foliorum, Moraxella_catarrhalis, Mycobacterium_smegmatis, Paenibacillus_larvae, Pectobacterium_atrosepticum, Propionibacterium_acnes, Proteus_mirabilis, Pseudomonas_aeruginosa, Pseudomonas_fluorescens, Pseudomonas_syringae, Ralstonia_solanacearum, Rhodococcus_erythropolis, Salmonella_enterica, Salmonella_typhimurium, Shigella_flexneri, Shigella_sonnei, Staphylococcus_aureus, Staphylococcus_epidermidis, Streptococcus_pneumoniae, Streptococcus_thermophilus, Streptomyces_griseus, Streptomyces_lividans, Streptomyces_venezuelae, Sulfolobus_islandicus, Vibrio_anguillarum, Vibrio_cholerae, Vibrio_parahaemolyticus, Yersinia_enterocolitica
```
Target genera
```
Achromobacter, Acidianus, Acinetobacter, Aeromonas, Arthrobacter, Bacillus, Brevibacillus, Brucella, Burkholderia, Campylobacter, Caulobacter, Cellulophaga, Citrobacter, Clostridium, Corynebacterium, Cronobacter, Cutibacterium, Dickeya, Enterobacter, Enterococcus, Erwinia, Escherichia, Flavobacterium, Gordonia, Haloarcula, Halorubrum, Helicobacter, Klebsiella, Lactobacillus, Lactococcus, Leuconostoc, Listeria, Mannheimia, Microbacterium, Moraxella, Mycobacterium, Paenibacillus, Pectobacterium, Prochlorococcus, Propionibacterium, Proteus, Pseudoalteromonas, Pseudomonas, Ralstonia, Rhizobium, Rhodococcus, Ruegeria, Salmonella, Serratia, Shigella, Sinorhizobium, Staphylococcus, Stenotrophomonas, Streptococcus, Streptomyces, Sulfolobus, Synechococcus, Thermus, Vibrio, Xanthomonas, Xylella, Yersinia
```


### Scripts
Main script:
   * **vHULK_v0.1.py** - Deep learning host prediction from phage genomes
  
Auxiliary script:
   * **download_and_set_models.py** - Used only at the first use for downloading and setting models
   
### Dependencies

All scripts from this project were coded in [Python 3](https://www.python.org/). So, first of all, make sure you have it installed and updated.  
v.HULK's main scrip (vHULK.py) requires Prokka and HMMER tools as dependencies. By installing Prokka and its dependencies, you will usually install HMMER tools automatically.  

* [Prokka](https://github.com/tseemann/prokka) - Rapid Prokaryotic genome annotation.
* [HMMER Tools](http://www.hmmer.org/) - Biosequence analysis using profile hidden Markov models

These Python libraries are required:

* [Numpy](http://www.numpy.org/), [Pandas](https://pandas.pydata.org/), [Scipy](https://www.scipy.org/) - Efficiently handling arrays and scientific computing for Python
* [Biopython](http://biopython.org/) - Handling biological sequences
* [Tensorflow](https://www.tensorflow.org/) - Google's Deep Neural Networks libraries 


To install these Python libraries, there are usually two ways: pip or conda.  
We strongly recomend creation of a specific conda environment for running v.HULK containing the installed dependencies.  
For example:

```
conda create --name vHULK python=3.6 numpy scipy biopython tensorflow pandas
conda activate vHULK
```

### Installing

Getting v.HULK ready to run is as simple as cloning this Github project or download and extract it to a directory inside your computer:

```
git clone https://github.com/LaboratorioBioinformatica/v.HULK
```

### Quick start

Inside the directory where v.HULK was extracted (or cloned), you will need to download and set the models. 
This is required only once and it is simple. Just run:
```
python3 download_and_set_models.py
```
All set!  
Now, to run v.HULK type:
```
python3 vHULK_v0.1.py -i input_directory -t num_threads
```

Change 'input_directory' to the folder where bins or genomes are stored in fasta format and 'num_threads' to the number of CPU cores to be used. Several threads should be used to speed up prokka and hmm searches.  
Results will be stored in the 'Results' folder inside the input directory.  
Obs: You need to execute the scripts from the directory where v.HULK was extracted, i.e., v.HULK's root folder. 

### Running the example datasets

We provide a folder with example datasets containing mocking bins of RefSeq viral and bacterial genomes.  
To try these examples, run:

```
python3 vHULK_v0.1.py -i test_input/ -t 12
```
