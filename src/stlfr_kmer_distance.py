#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 08:06:52 2019

@author: Zewei Song
@email: songzewei@genomics.cn
"""
from __future__ import print_function
from __future__ import division
import argparse
import json
from metaSeq import bead
from metaSeq import kmer
from itertools import combinations

parser = argparse.ArgumentParser()
parser.add_argument('-i', help='Input JSON-Bead file.')
parser.add_argument('-t', '--threshold', nargs='+', default=[0.02,0.04], type = float, help='Threshold for distance.')
parser.add_argument('-rawout', help='Output raw edge file.')
parser.add_argument('-tout', help='Output thresholded edge file.')
parser.add_argument('-k', default=21, type=int, help='Kmer size default = 21')
args = parser.parse_args()

inputFile = args.i
outputRaw = args.rawout
outputThreshold = args.tout
threshold = args.threshold
kmerSize = args.k

# Read in JSON-Bead file
#Calculate kmer pools for all beads
kmerPool = []
beadCount = 0
with open(inputFile, 'r') as f:
    for line in f:
        b = bead.beadSequence(json.loads(line))
        kmerPool.append(kmer.kmerCount(b, kmerSize))
        beadCount += 1
print('Found {0} beads.'.format(beadCount))

# Calculate kmer distance for all pairs
edge = []
edgeThreshold = []
n1 = 0
n2 = 0
for pair in combinations(kmerPool, 2):
    D = kmer.kmerDistance((pair[0].set, pair[1].set)).mashDistance()
    edge.append((pair[0].barcode, pair[1].barcode, D))
    if threshold[0] <= D <= threshold[1]:
        edgeThreshold.append((pair[0].barcode, pair[1].barcode, D))
        n2 += 1
    n1 += 1
    print('Calculated {0} pairs, {1} fall into threshold.'.format(n1, n2))

# Write to output
if outputRaw:
    with open(outputRaw, 'w') as f:
        f.write('Source\tTArget\tDistance\n')
        for line in edge:
            f.write('{0}\t{1}\t{2}\n'.format(line[0], line[1], line[2]))
with open(outputThreshold, 'w') as f:
    f.write('Source\tTArget\tDistance\n')
    for line in edgeThreshold:
        f.write('{0}\t{1}\t{2}\n'.format(line[0], line[1], line[2]))
