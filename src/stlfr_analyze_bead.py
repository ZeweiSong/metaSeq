#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 16:47:24 2018

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

# Read in the JSON bead file
inputFile = 'mock.qc.derep.json'
beadList = []
with open(inputFile, 'r') as f:
    for line in f:
        beadList.append(bead.beadSequence(json.loads(line)))
#%%
# Write the first 10 bead into FASTA file
i = 1
threshold = 100
for line in beadList:
    #outputFile = line.barcode + '.fa'
    line.fastaWrite()
    i += 1
    if i > 10:
        break