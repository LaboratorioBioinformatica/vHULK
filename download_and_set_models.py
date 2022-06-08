#!/usr/bin/python3
# coding: utf-8
#
#
## This is v.HULK, a toolkit to find phage hosts
## This is an auxiliary script to download and set models
## Developed by Deyvid Amgarten

# Libraries
import os
import subprocess

# Greeting message
print('\nYou only need to run this script once!\n')

# Verify and download DBs
if not os.path.isdir('models'):
    print('Downloading trained models. Do not worry, that will just take a few minutes and is executed only in the first time... \n')
    os.system('wget http://projetos.lbi.iq.usp.br/phaghost/vHULK/models/models_Jun_2022.tar.xz')
    print('Extracting models files...\n')
    if subprocess.call('tar -xzf models_Jun_2022.tar.xz', shell=True) == 1:
        print('Error extracting models. Please this script again\n')
        quit()
    print('Trained models are all set!\n')
else:
    print('Trained models are already set.\n')

if not os.path.isfile('models/all_vogs_hmm_profiles_feb2018.hmm.h3m'):
    print('Downloading flat file database. Do not worry, that will just take a few minutes and is executed only in the first time... \n')
    os.system('wget http://projetos.lbi.iq.usp.br/phaghost/vHULK/models/AllvogHMMprofiles.tar.gz')
    print('Extracting database file...\n')
    if subprocess.call('tar -xzf AllvogHMMprofiles.tar.gz', shell=True) == 1:
        print('Error extracting database\n')
        quit()
    subprocess.run('cat AllvogHMMprofiles/* > models/all_vogs_hmm_profiles_feb2018.hmm', shell=True)
    subprocess.run('rm -r AllvogHMMprofiles/ AllvogHMMprofiles.tar.gz', shell=True)
    print('Compressing hmm database...')
    if subprocess.call('hmmpress models/all_vogs_hmm_profiles_feb2018.hmm', shell=True) == 1:
        print('Error using hmmer tools (hmmpress). Verify if it is installed!\n')
        quit()
    print('Database is all set!\n')
else:
    print('HMM Database is already set.\n')

if not os.path.isdir('datasets'):
    print('Downloading datasets. Do not worry, that will just take a few minutes and is executed only in the first time... \n')
    os.system('wget http://projetos.lbi.iq.usp.br/phaghost/vHULK/datasets/training_dataset_Jun_2022.tar.xz')
    print('Extracting models files...\n')
    if subprocess.call('tar -xzf training_dataset_Jun_2022.tar.xz', shell=True) == 1:
        print('Error extracting models. Please this script again\n')
        quit()
    print('Trained models are all set!\n')
else:
    print('Trained models are already set.\n')        
    
print('Thank you for using v.HULK.')
