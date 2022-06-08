import pandas as pd
import numpy as np
import csv
import re

#define here the path for aai matrix
df_genus = pd.read_csv('aai_summary.tsv', sep='\t')
dic_rep = {}
dic_exc = {}
count = 0

# Set here the threshold of AAI to assembly clusters and collapse sequenctially and also take into account the threshold of
# orthologous fraction between genomes (>70%)
for index, row in df_genus.loc[(df_genus['Orthologous fraction (OF)'] >= 70) & (df_genus['Mean AAI'] >= 70)].iterrows():
    count += 1
    if row['#Genome A']+'|'+str(row['Genes in A']) not in dic_rep:
        flag = 0
        for kn in dic_rep:
            for knl in dic_rep[kn]:
                if row['#Genome A']+'|'+str(row['Genes in A']) == knl:
                    flag = 1
                    break
        if flag == 0:
            dic_rep[row['#Genome A']+'|'+str(row['Genes in A'])] = [row['Genome B']+'|'+str(row['Genes in B'])]
    else:
        dic_rep[row['#Genome A']+'|'+str(row['Genes in A'])].append(row['Genome B']+'|'+str(row['Genes in B']))
    if count%100 == 0:
        print(count) 


#Verify clusters, choose the longest representative and take host into account, leaving one representative per host
list_rep = []
list_exc = []
for key in dic_rep:
    # Number of genes
    maxi = int(key.split('|')[-1])
    # Format: ACC-host|#genes
    maxi_host = key.split('|')[0].split('-')[-1]
    maxi_key = key
    # List of unique hosts
    host_control = [key]
    # Iterate list of genomes similar to the key
    for item in dic_rep[key]:
        # Verify if #genes do genoma da lista é maior do que o da key
        if item.split('|')[0].split('-')[-1] == maxi_host:
            
            if int(item.split('|')[-1]) > maxi:
                # In this case, exclude key
                list_exc.append(maxi_key)
                # Define new values to maxi
                maxi = int(item.split('|')[-1])
                maxi_key = item
            else:
                # in this case, exclude genome from list
                list_exc.append(item)
        else:
            if len(host_control) == 1:
                host_control.append(item)
            else:
                flagi = 0
                for host in host_control:
                    if item.split('|')[0].split('-')[-1] == host.split('|')[0].split('-')[-1]:
                        if  int(item.split('|')[-1]) < int(host.split('|')[-1]):
                            list_exc.append(item)
                        else:
                            list_exc.append(host)
                        flagi == 1
                if flagi == 0:
                    host_control.append(item)
                                       
    list_rep.append(maxi_key)


# Exclude redundant genomes from the baseline dataset
for each in set(list_exc):
    #print (each)
    file = each.split('|')[0]+'.fasta'
    print(file)
    # Path to baseline FASTA genomes
    !rm redundancy_removed/genus/fasta_genus_AAI70/{file}