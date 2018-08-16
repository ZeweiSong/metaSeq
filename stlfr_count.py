#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 15 17:37:43 2018

@author: Zewei Song
@email: songzewei@genomics.cn
"""
from __future__ import print_function
from __future__ import division
from metaSeq import io as seqIO
from metaSeq import bead
import json
import random
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', help='Input beadJson file')
parser.add_argument('-o', help='Output short fragment count per bead')
parser.add_argument('-d', help='Output distribution of count')
args = parser.parse_args()

inputFile = args.i
outputFile = args.o
outputDist = args.d

count = {}
dist = {}
with open(inputFile, 'r') as f:
    for line in f:
        currentBead = bead.beadSequence(json.loads(line.strip('\n')))
        fragmentCount = len(currentBead.fragments)
        count[currentBead.barcode] = fragmentCount
        try:
            dist[fragmentCount] += 1
        except KeyError:
            dist[fragmentCount] = 1

# Output
with open(outputDist, 'w') as f:
    f.write('{0}\t{1}\n'.format('FragmentCount','Frequency'))
    for key, value in dist.items():
        f.write('{0}\t{1}\n'.format(str(key), str(value)))

with open(outputFile, 'w') as f:
    f.write('{0}\t{1}\n'.format('Barcode','FragmentCount'))
    for key, value in count.items():
        f.write('{0}\t{1}\n'.format(key, str(value)))