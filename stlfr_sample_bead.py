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
args = parser.parse_args()

inputFile = args.i.split(',')
outputFile = args.o
n = args.n

beadJson = seqIO.fastq2json(inputFile[0], inputFile[1], barcodePosition=-1)
print(len(beadJson))
# print(beadJson[list(beadJson.keys())[0]])
# Get the n barcode randomly
randomKeys = random.sample(list(beadJson.keys()), n)

with open(outputFile, 'w') as f:
    for key in randomKeys:
        f.write('{0}\n'.format(json.dumps({key:beadJson[key]})))