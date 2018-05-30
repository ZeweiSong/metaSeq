#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 16:59:32 2017

Parse the NCBI genome summary, and UNITE taxon

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%%
from __future__ import print_function
from __future__ import division
from metaSeq import io


# Read in NCBI and UNITE list
# Seaerch UNITE for the full taxonomy of NCBI record
#path = 'E:/LabSpace/ncbi_parser'
#os.chdir(path)
#print('Set working folder to {0}'.format(path))

ncbi_file = 'fungi_assembly_summary.txt' # Download from NCBI
unite_file = 'sh_general_release_dynamic_01.12.2017.fasta' # UNITE

ncbi_fungi = {}
unite_species = {} # Dictionary for taxonomy using species as index
unite_genus = {} # Dictionary for taxonomy using Genus as index
print('Use {0} as UNITE database.'.format(unite_file))
unite = io.sequence(unite_file)
for record in unite:
    label = record[0].split('|')
    genus = label[0].split('_')[0]
    unite_species[label[0]] = label[-1]
    unite_genus[genus] = ';'.join(label[-1].split(';')[:6])
print('\tFound {0} fungal Species in UNITE.'.format(len(unite_species)))
print('\tFound {0} fungal Genus in UNITE.'.format(len(unite_genus)))

print('Read in NCBI Fungl Ref Genome list: {0}'.format(ncbi_file))
with open(ncbi_file, 'r') as f:
    comment = f.readline()
    title = f.readline()
    title = title[2:].strip('\n').split('\t')
    for line in f:
        line = line.strip('\n').split('\t')
        ncbi_fungi[line[0]] = {}
        for key, item in zip(title, line):
            ncbi_fungi[line[0]][key] = item
print('\tFound {0} record in NCBI Fungal Ref Genome list.'.format(len(ncbi_fungi)))

not_found = []
species_count = 0
genus_count = 0
for key, value in ncbi_fungi.items():
    species_name = value['organism_name'].split(' ')[1]
    genus = value['organism_name'].split(' ')[0]
    species = genus + '_' + species_name
    try:
        ncbi_fungi[key]['taxonomy'] = unite_species[species]
        species_count += 1
    except KeyError: # No current species in UNITE, try to find the same Genus
        try:
            ncbi_fungi[key]['taxonomy'] = unite_genus[genus] + ';' + species
            genus_count += 1
        except KeyError:
            not_found.append([key, species])
print('\t\t{0} records have corresponding Species in UNITE.'.format(species_count))
print('\t\t{0} records have corresponding Genus in UNITE.'.format(genus_count))
print('\t\t{0} record cannot be found in UNITE.'.format(len(not_found)))
#%% Read in possible file suffix
suffix = []
suffix_file = 'assembly_file_list.txt'
with open(suffix_file, 'r') as f:
    for line in f:
        suffix.append(line.strip('\n'))
print('Genome sequence: {0}'.format(suffix[5]))
print('Genome assembly report: {0}'.format(suffix[0]))
print('Genome assembly statistics: {0}'.format(suffix[1]))
#%% Test on several genomes
import os
from urllib import request
names = list(ncbi_fungi.keys())
print('Find {0} genomes'.format(len(names)))
suffix_index = [0,1,5] # Download files of # 0, 1, 4
# Create new folder for each genome
for genome in names:
    if not os.path.exists(genome):
        os.mkdir(genome)

for genome in names:
    ftp_path = ncbi_fungi[genome]['ftp_path']
    prefix = ncbi_fungi[genome]['ftp_path'].split('/')[-1]
    saved_files = []
    download_links = []
    for i in suffix_index:
        saved_files.append(genome + '/' + prefix + '_' + suffix[i])
        download_links.append(ftp_path + '/' + prefix + '_' + suffix[i])
    #print(saved_files)
    #print(download_links)
#%%
import time
for file, link in zip(saved_files, download_links):
    t1 = time.time()
    print('Downloading from: {0} ...'.format(link))   
    request.urlretrieve(link, file)
    t2 = time.time()
    print('Saved to {0}, used {1} seconds.'.format(file, t2 - t1))