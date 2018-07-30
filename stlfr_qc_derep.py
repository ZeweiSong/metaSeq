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

inputFile = 'mock.json'
maxee = 1
outputFile = 'mock.qc.derep.json'

mockRaw = []
with open(inputFile, 'r') as f:
    for line in f:
        mockRaw.append(json.loads(line))

#%% For each bead, remove low quality read and duplicated reads
# Then write to a new JSON file
for item in mockRaw:
    beadQC = bead.maxEE(item)
    beadDerep = bead.derep(beadQC)
    beadProcessed = bead.beadSequence(beadDerep)
    if len(beadProcessed.fragments) > 0:
        beadProcessed.jsonWrite(outputFile, mode='a')

#%% The JSON file can be read in by line
# A single line can be converted to a bead Class
beadList = []
with open(outputFile, 'r') as f:
    for line in outputFile:
        beadList.append(bead.beadSequence(json.loads(line)))