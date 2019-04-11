#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 16:56:16 2018

@author: Zewei Song
@email: songzewei@genomics.cn
"""
from __future__ import print_function
from __future__ import division
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='Input alignments')
parser.add_argument('-o', '--output', help='Output paired alignments')
args = parser.parse_args()
inputAln = args.input.split(',')
outputAln = args.output

# Build the query dictionary
# Two alignments are combined into one
#%%
#inputAln = ['501.r1fw.b6', '501.r2rv.b6']
queryDict = {}
with open(inputAln[0], 'r') as f:
    for line in f:
        line = line.strip('\n').split('\t')
        line[0] = line[0].split('/')[0]
        queryDict[line[0]] = {'r1':[], 'r2':[]}
        queryDict[line[0]]['r1'].append(line[1])

with open(inputAln[1], 'r') as f:
    for line in f:
        line = line.strip('\n').split('\t')
        line[0] = line[0].split('/')[0]
        try:
            queryDict[line[0]]['r2'].append(line[1])
        except KeyError:
            queryDict[line[0]] = {'r1':[], 'r2':[]}
            queryDict[line[0]]['r2'].append(line[1])
#%%
aln = {}
for query, alignments in queryDict.items():
    try:
        overlap = tuple(set(alignments['r1']).intersection(set(alignments['r2'])))
        if len(overlap) > 0:
            aln[query] = overlap
    except KeyError:
        pass
#%%
with open(outputAln, 'w') as f:
    for query, alignments in aln.items():
        for item in alignments:
            f.write('{0}\t{1}\n'.format(query, item))