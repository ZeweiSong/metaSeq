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
queryDict = {}
with open(inputAln[0], 'r') as f:
    for line in f:
        line = line.strip('\n').split('\t')
        try:
            queryDict[line[0]]['r1'].append(line[1])
        except KeyError:
            try:
                queryDict[line[0]]['r1'] = [line[1]]
            except KeyError:
                queryDict[line[0]] = {}
with open(inputAln[1], 'r') as f:
    for line in f:
        line = line.strip('\n').split('\t')
        try:
            queryDict[line[0]]['r2'].append(line[1])
        except KeyError:
            try:
                queryDict[line[0]]['r2'] = [line[1]]
            except KeyError:
                queryDict[line[0]] = {}
aln = {}
for query, alignments in queryDict.items():
    overlap = tuple(set(alignments[0]).intersection(set(alignments[1])))
    if len(overlap) > 0:
        aln[query] = overlap

with open(outputAln, 'w') as f:
    for query, alignments in aln.items():
        for item in alignments:
            f.write('{0}\t{1}\n'.format(query, item))