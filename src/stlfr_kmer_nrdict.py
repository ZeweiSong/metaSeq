#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 10:36:47 2018

Create a non-redundant kmer dictionary linked to bead barcodes.

Using dictionary turns out to be a faster way to calculate pairwise distance
than list overlap methods.

This script is under development.

Coders who love to comment their code are unlikely to have bad luck.

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%%
from __future__ import print_function
from __future__ import division
import argparse
import json
from metaSeq import io as seqIO

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='Input Kmer json file.')
parser.add_argument('-o', '--output', help='Output file.')
args = parser.parse_args()
args = parser.parse_args(['-i', 'kmer.sample.json', '-o', 'kmer.sample.nrdt.json'])
inputFile = args.input
outputFile = args.output

nrd = {}
count = 0
for item in seqIO.beadJson(inputFile):
    count += 1
    barcode = list(item.keys())[0]
    kmers = item[barcode]
    for km in kmers:
        nrd[km] = nrd.get(km, {})
        nrd[km].update({barcode:1})
    print(count, len(nrd))
    
print('Found {0} non-redundant kmers.'.format(len(nrd)))
with open(outputFile, 'w') as f:
    f.write('{0}\n'.format(json.dumps(nrd)))


