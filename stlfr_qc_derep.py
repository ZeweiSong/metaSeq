#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 15:42:25 2018

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%%
from __future__ import print_function
from __future__ import division
import json
from metaSeq import io as seqIO
from metaSeq import qc as seqQC
from metaSeq import bead
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', help='Input JSON bead file.')
parser.add_argument('-maxee', default=1, type=int, help='Threshold of maxEE.')
parser.add_argument('-o', help='Output JSON bead file.')
args = parser.parse_args()

inputFile = args.i
maxee = args.maxee
outputFile = args.o

mockRaw = []
with open(inputFile, 'r') as f:
    for line in f:
        mockRaw.append(json.loads(line))

#%% For each bead, remove low quality read and duplicated reads
# Then write to a new JSON file
for item in mockRaw:
    beadQC = bead.maxEE(item, maxee=maxee)
    beadDerep = bead.derep(beadQC)
    beadProcessed = bead.beadSequence(beadDerep)
    if len(beadProcessed.fragments) > 0:
        beadProcessed.jsonWrite(outputFile, mode='a')
#%% The JSON file can be read in by line
# A single line can be converted to a bead Class
'''
beadList = []
with open(outputFile, 'r') as f:
    for line in f:
        beadList.append(bead.beadSequence(json.loads(line)))
'''