#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 16:37:47 2019

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%%
from __future__ import print_function
from __future__ import division
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='Input alignment')
parser.add_argument('-o', '--output', help='Output report file')
args = parser.parse_args()
args = argparse.Namespace(input ='test.b6', output = 'test.txt' )
inputAln = args.input
outputReport = args.output

alnDict = {}
refDict = {}
with open(inputAln, 'r') as f:
    for line in f:
        line = line.strip('\n').split('\t')
        query = line[0]
        ref = line[1]
        try:
            alnDict[query] += 1
        except KeyError:
            alnDict[query] = 1
        try:
            refDict[ref] += 1
        except KeyError:
            refDict[ref] = 1

queryDist = {}
refDist = {}

for key, queryAlnCount in alnDict.items():
    try:
        queryDist[queryAlnCount] += 1
    except KeyError:
        queryDist[queryAlnCount] = 1
maxAlnCount = max([i for i in queryDist.keys()])
qdist = []
for i in range(maxAlnCount):
    try:
        qdist.append((i+1, queryDist[i+1]))
    except KeyError:
        qdist.append((i+1, 0))

for key, refAlnCount in refDict.items():
    try:
        refDist[refAlnCount] += 1
    except KeyError:
        refDist[refAlnCount] = 1
maxAlnCount = max([i for i in refDist.keys()])
rdist = []
for i in range(maxAlnCount):
    try:
        rdist.append((i+1, refDist[i+1]))
    except KeyError:
        rdist.append((i+1, 0))

with open(outputReport, 'w') as f:
    f.write('QueryAlnCount\tCount\n')
    for line in qdist:
        f.write('{0}\t{1}\n'.format(line[0], line[1]))
    f.write('RefAlnCount\tCount\n')
    for line in rdist:
        f.write('{0}\t{1}\n'.format(line[0], line[1]))