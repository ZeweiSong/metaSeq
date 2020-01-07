#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 10:36:47 2018

Coders who love to comment their code are unlikely to have bad luck.

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%%
from __future__ import print_function
from __future__ import division
import json
inputBiom = 'its2.taxa.biom'
with open(inputBiom, 'r') as f:
    table = json.load(f)
print('Read in {0} as biom-format JSON version 1.0.'.format(inputBiom))
print('Format: {0}'.format(table['format']))
print('Matrix type: {0}'.format(table['matrix_type']))
print('Matrix has {0} rows/features, and {1} columns/samples.'.format(table['shape'][0], table['shape'][1]))

# Build a new format for the table
# We know you hate biom format, so we bring you zeta-format!
# Zeta-Format is a strict JSON format, aka, just read'em in as Python dictionary.
# For a new table, refer species using 'bugs'
# refer samples using 'samples'
# All bugs are keyed using their index in the biom file,
# This index is corresponding to the abundance dict in samples, so number, instead of species ids, are used repeatedly.
new_table = {'bugs':{}, 'samples':{}}
# Add rows, use i as the key
for i, item in enumerate(table['rows']):
    new_table['bugs'][i] = item
    new_table['bugs'][i]['metadata']['taxonomy'] = ';'.join(new_table['bugs'][i]['metadata']['taxonomy'])
# Add samples, use sample id as the key
for item in table['columns']:
    new_table['samples'][item['id']] = item
    new_table['samples'][item['id']]['abundance'] = {}
# Add abundance data to samples
for item in table['data']:
    sample_id = table['columns'][item[1]]['id']
    new_table['samples'][sample_id]['abundance'][item[0]] = item[2]

#%% We can now slice the table using samples id!
sample_list = tuple(new_table['samples'].keys())
sample_yc = [i for i in sample_list if i[0:2] in ['XX','YC','WS','SZ']]
table_yc = {'bugs':{}, 'samples':{}}
bug_yc = []
for item in sample_yc:
    table_yc['samples'][item] = new_table['samples'][item]
    bug_yc += list(table_yc['samples'][item]['abundance'].keys())
bug_yc = tuple(set(bug_yc))
for item in bug_yc:
    table_yc['bugs'][item] = new_table['bugs'][item]

#%% Now we would like to write it as a tsv table
bug_list = tuple(table_yc['bugs'].keys())
sample_list = tuple(table_yc['samples'].keys())
table_header = ['BUGS_ID'] + list(table_yc['samples'].keys()) + list(table_yc['bugs'][bug_list[0]]['metadata'].keys())
table_tsv = []
for bug in bug_list:
    bug_id = table_yc['bugs'][bug]['id']
    line = [bug_id]
    for sample in sample_list:
        abundance_value = table_yc['samples'][sample]['abundance'].get(bug, 0)
        line.append(abundance_value)
    for key, value in table_yc['bugs'][bug]['metadata'].items():
        line.append(value)
    table_tsv.append(line)
outputTsv = 'caas_its2_2018.tsv'
table_tsv.sort(key=lambda x:sum(x[1:-1]), reverse=True)
with open(outputTsv, 'w') as f:
    f.write('{0}\n'.format('\t'.join(table_header)))
    for item in table_tsv:        
        f.write('{0}\n'.format('\t'.join([str(i) for i in item])))

