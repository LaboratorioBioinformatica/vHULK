# v.HULK
## Viral Host UnveiLing Kit  
A tool-kit for phage host prediction


### Scripts
Main script:
   * **vHULK_v0.1.py** - Deep learning host prediction from phage genomes
  
Auxiliary script:
   * **download_and_set_models.py** - Used only at the first use for downloading and setting models
   
## Dependencies

All scripts from this project were coded in [Python 3](https://www.python.org/). So, first of all, make sure you have it installed and updated.  
v.HULK's main scrip (vHULK.py) requires Prokka and HMMER tools as dependencies. By installing Prokka and its dependencies, you will usually install HMMER tools automatically.  

* [Prokka](https://github.com/tseemann/prokka) - Rapid Prokaryotic genome annotation.
* [HMMER Tools](http://www.hmmer.org/) - Biosequence analysis using profile hidden Markov models

These Python libraries are required:

* [Numpy](http://www.numpy.org/), [Scipy](https://www.scipy.org/), [Pandas]() - Efficiently handling arrays and scientific computing for Python
* [Biopython](http://biopython.org/) - Handling biological sequences
* [Tensorflow](https://www.tensorflow.org/) - Google's Deep Neural Networks libraries 


To install these Python libraries, there are usually two ways: pip or conda.  
We strongly recomend creation of a specific conda environment for running v.HULK containing the installed dependencies.  
For example:

```
conda create --name vHULK python=3.6 numpy scipy biopython tensorflow pandas
```

## Installing

Getting v.HULK ready to run is as simple as cloning this Github project or download and extract it to a directory inside your computer:

```
git clone https://github.com/LaboratorioBioinformatica/v.HULK
```

## Quick start

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

## Running the example datasets

We provide a folder with example datasets containing mocking bins of RefSeq viral and bacterial genomes.  
To try these examples, run:

```
python3 vHULK_v0.1.py -i example_data/bins_8k_refseq -t 12
```
