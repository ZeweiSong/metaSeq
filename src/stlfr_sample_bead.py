#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 15 14:25:43 2018

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%%
from __future__ import print_function
from __future__ import division
from metaSeq import io as seqIO
import json
import random
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', help='Input file to sample')
parser.add_argument('-o', help='Output file')
parser.add_argument('-n', type=int, help='Number of beads to sample')
parser.add_argument('-not_gz', action='store_false', help='Specify if the input is not a gz file, in most case you do not need it.')
args = parser.parse_args()

inputFile = args.i
outputFile = args.o
n = args.n
not_gz = args.not_gz

beadJson = []
with open(inputFile, 'r') as f:
    for line in f:
        beadJson.append(line)

print('Found {0} beads.'.format(len(beadJson)))
# print(beadJson[list(beadJson.keys())[0]])
# Get the n barcode randomly
randomBeads = random.sample(beadJson, n)

with open(outputFile, 'w') as f:
    for line in randomBeads:
        f.write('{0}'.format(line))
print('Randomly sampled {0} beads into {1}.'.format(n, outputFile))