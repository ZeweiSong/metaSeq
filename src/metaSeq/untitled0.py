#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 10:36:47 2018

Parse the GenBank format into orgnized FASTA database.

Coders who love to comment their code are unlikely to have bad luck.

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%%
from __future__ import print_function
from __future__ import division
from Bio import SeqIO
from metaSeq import io as metaSeqIO

def parse_level(l, d):
    levels = ('genus', 'family', 'order', 'class', 'phylum', 'superkingdom')
    taxa = {'genus':'Unassigned', 'family':'Unassigned', 'order':'Unassigned', \
            'class':'Unassigned', 'phylum':'Unassigned', 'superkingdom':'Unassigned'}
    for item in l:
        value = d.get(item, 'Unknown')
        if value in levels:
            taxa[value] = item
    return taxa
 
#%% 2019/9/11 Read in nodes.dmp, and use the taxid as key, and taxa level as value
countLine = 0
levelSpace = {}
taxid_level = {}
with open('E:/Bench/taxdump/nodes.dmp', 'r') as f:
    for line in f:
        countLine += 1
        line = line.strip('\n').split('\t|\t')
        child = int(line[0])
        parent = int(line[1])
        levelSpace[line[2]] = levelSpace.get(line[2], 0) + 1
        taxid_level[child] = line[2]
print('nodes.dmp has {0} lines.'.format(countLine))
print(levelSpace)

#% Read in names.dmp
names = {0:[('Unclassified', '', 'scientific name')]}
with open('E:/Bench/taxdump/names.dmp', 'r') as f:
    for line in f:
        line = line.strip('\t|\n').split('\t|\t')
        nodeid = int(line[0])
        value = line[1]
        unique_name = line[2]
        name_class = line[3]
        names[nodeid] = names.get(nodeid, []) + [(value, unique_name, name_class)]
print('Found {0} names in names.dmp.'.format(len(names)))

# Clean for scientific names
scientific_names = {}
for key, value in names.items():
    for record in value:
        if record[2] == 'scientific name':
            scientific_names[key] = scientific_names.get(key, []) + [record]

#% 2019/9/11 Use scientific name string as key, and level as value
string_level = {}
for key, value in scientific_names.items():
    if key != 0:
        string_level[value[0][0]] = taxid_level[key]
string_level['Methanosaetaceae'] = 'family'
print(len(string_level))
#%%
gb_file = '33175_33317.gb'
gb_record = []
count = 0
tax = []
for record in SeqIO.parse(open(gb_file, 'r'), 'genbank'):
    seq = record.seq.__str__()
    annot = record.annotations
    organism = '_'.join(annot['organism'].split(' '))
    species = organism.split('_')[1]
    taxa = parse_level(annot['taxonomy'], string_level)
    taxonomy = 'k__{0};p__{1};c__{2};o__{3};f__{4};g__{5};s__{6}'.format(taxa['superkingdom'], taxa['phylum'], taxa['class'], \
                   taxa['order'], taxa['family'], taxa['genus'], taxa['genus'] + '_' + species)
    head = annot['accessions'][0] + '.' + str(annot['sequence_version']) + '|' + organism + '|' + taxonomy
    gb_record.append((head, seq))
    count += 1
print(count)

#%% Write
count = metaSeqIO.write_seqs(gb_record, 'ncbi_targeted_refseq_16s.fa', fastx='a', gz=False)
print(count)
#%% Also write to the format of Burst
part1 = []
part2 = []
for record in gb_record:
    head = record[0].split('|')
    part1.append('{0}\t{1}\n'.format(head[0] + '|' + head[1], head[2]))
    part2.append((head[0] + '|' + head[1], record[1]))
count = metaSeqIO.write_seqs(part2, 'ncbi_targeted_refseq_16s_burst.fa', fastx='a', gz=False)
with open('ncbi_targeted_refseq_16s_burst.txt', 'w') as f:
    for item in part1:
        f.write(item)

#%% Check for names not in names.dmp and nodes.dmp, there are some
c = {}
unassign = []
for item in tax:
    for value in item:
        try:
            c[string_level[value]] = c.get(string_level[value], 0) + 1
        except KeyError:
            unassign.append(value)
            print(item)