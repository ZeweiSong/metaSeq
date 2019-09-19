#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 10:36:47 2018

Calculate Kmer distance from kmer json file

Coders who love to comment their code are unlikely to have bad luck.

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%%
from __future__ import print_function
from __future__ import division
import argparse
from metaSeq import kmer
from metaSeq import io as seqIO
from itertools import combinations
import random

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='JSON kmer file.')
parser.add_argument('-o', '--output', help='Output file.')
args = parser.parse_args()
# args = parser.parse_args(['-i', 'kmer.sample.json', '-o', 'kmer.distance.tsv'])
inputFile = args.input
outputFile = args.output

print('Read in kmer json file')
kmerParser = seqIO.beadJson('kmer.sample.json')
kmers = []
for item in kmerParser:
    currentBead = list(list(item.items())[0])
    currentBead[1] = random.sample(currentBead[1], len(currentBead[1])//10) 
    kmers.append(currentBead)
print('Found {0} beads'.format(len(kmers)))
print('Start calculating kmer distance')
mashD = []
count = 0
for pairs in combinations(kmers, 2):
    count += 1
    if count // 100000 > 0 and count % 10000 == 0:
        print(count)
    k1 = pairs[0]
    k2 = pairs[1]
    kd = kmer.kmerDistance((k1[1], k2[1]))
    mashD.append((k1[0], k2[0], kd.mashDistance()))
with open(outputFile, 'w') as f:
    f.write('{0}\t{1}\t{2}\n'.format('Target', 'Source', 'mashDistance'))
    for line in mashD:
        f.write('{0}\t{1}\t{2}\n'.format(line[0],line[1],line[2]))