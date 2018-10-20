#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 09:59:28 2018

@author: Zewei Song
@email: songzewei@genomics.cn
"""
from __future__ import print_function
from __future__ import division
from metaSeq import io as seqIO
from metaSeq import bead
from metaSeq import kmer
import json
from itertools import combinations

#%% Try to sample bead with more than 1 fragments
inputBeadJson = 'even10pg_10over4.cutmerge.qcderep.json'

beadPool = []
with open(inputBeadJson, 'r') as f:
    for line in f:
        beadPool.append(bead.beadSequence(json.loads(line)))
print('Find {0} beads.'.format(len(beadPool)))

#%% Get beads with more than 1 fragment
beadPool_1 = []
for item in beadPool:
    if len(item.fragments) > 1:
        beadPool_1.append(item)
print('Find {0} bead with more than 1 fragments.'.format(len(beadPool_1)))

#%% Write these bead to FASTA files
outputFolder = 'beadPool_1/'

for item in beadPool_1:
    item.fastaWrite(folder=outputFolder)
