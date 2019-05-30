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
#args = argparse.Namespace(input ='A_R1_522_AB_FORAGE.b6,A_R2_522_AB_FORAGE.b6', output = 'test.b6' )
inputAln = args.input.split(',')
outputAln = args.output

# Build the query dictionary
# Two alignments are combined into one
#%%
#inputAln = ['501.r1fw.b6', '501.r2rv.b6']
queryDict = {}
print('Parsing Read1 alignment ...')
with open(inputAln[0], 'r') as f:
    for line in f:
        line = line.strip('\n').split('\t')
        line[0] = line[0].split('/')[0] # Get the Query label (stripped Read info)
        try:
            queryDict[line[0]]['r1'].append(line[1])
        except KeyError:
            queryDict[line[0]] = {'r1':[], 'r2':[]}
            queryDict[line[0]]['r1'].append(line[1])
print('Found {0} queries in Read1 alignment.'.format(len(queryDict)))
print('Parsing Read2 alignment ...')

r2inr1 = 0
r2notinr1 = 0
with open(inputAln[1], 'r') as f:
    for line in f:
        line = line.strip('\n').split('\t')
        line[0] = line[0].split('/')[0]
        try:
            queryDict[line[0]]['r2'].append(line[1])
        except KeyError:
            queryDict[line[0]] = {'r1':[], 'r2':[]}
            queryDict[line[0]]['r2'].append(line[1])
            r2notinr1 += 1
print('Total queries add to {0}'.format(len(queryDict)))
#%%
aln = {}
for query, alignments in queryDict.items():
    try:
        overlap = tuple(set(alignments['r1']).intersection(set(alignments['r2'])))
        if len(overlap) > 0:
            aln[query] = overlap
    except KeyError:
        pass

print('Found {0} queries in both alignment.'.format(len(aln)))
print('\t representing {0} alignments.'.format(sum([len(i) for i in aln.values()])))
print('Found {0} queries ONLY in Read2 alignment.'.format(r2notinr1))
#%%
with open(outputAln, 'w') as f:
    for query, alignments in aln.items():
        for item in alignments:
            f.write('{0}\t{1}\n'.format(query, item))