#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 10:36:47 2018

Coders who love to comment their code are unlikely to have bad luck.

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%
from __future__ import print_function
from __future__ import division
from __future__ import print_function
from __future__ import division
import argparse
import json
from metaSeq import io as seqIO

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='Input Kmer json file.')
parser.add_argument('-o', '--output', help='Output pairwise jacarrd.')
args = parser.parse_args()
#args = parser.parse_args(['-i', 'kmer.sample.json', '-o', 'kmer.sample.nrdt.json'])
inputFile = args.input
outputFile = args.output
with open('kmer.sample.nrdt.json', 'r') as f:
    nrd = json.load(f)

#% Test the pairwise distance
from itertools import combinations
beadPool = []
for item in seqIO.beadJson(inputFile):
    barcode = list(item.keys())[0]
    kmers = item[barcode]
    beadPool.append((barcode, kmers))

pd = []
count = 0
for pair in combinations(beadPool, 2):
    count += 1
    if len(pair[0]) < len(pair[1]):
        k1 = pair[0]
        k2 = pair[1]
    else:
        k1 = pair[1]
        k2 = pair[0]
    share = 0
    for item in k1[1]:
        share += nrd[item].get(k2[0], 0)
    jcd = share / (len(k1[1]) + len(k2[1]) - share)
    pd.append((k1[0], k2[0], jcd))
    print(count, count/29272726 * 100)
with open(outputFile, 'w') as f:
    f.write('{0}]\t{1}\t{2}\n'.format('Target', 'Source', 'Jacarrd'))
    for line in pd:
        f.write('{0}]\t{1}\t{2}\n'.format(line[0], line[1], line[2]))