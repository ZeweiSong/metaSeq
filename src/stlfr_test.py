#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 10:36:47 2018

Coders who love to comment their code are unlikely to have bad luck.

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%
from __future__ import print_function
from __future__ import division
from metaSeq import io as seqIO
from metaSeq import bead
from metaSeq import kmer
import json
import random
beads = seqIO.beadJsonTrunk('CL100077200_L01.json', trunk = 1000000)
report = []
k = 21
count = 0
# kmerDist = {}

print('Start counting kmers')
for trunk in beads:
    with open('kmerDist/dist_{0}.tsv'.format(count), 'w') as f:
        kmerDist = {}
        count += 1
        print(count)
        for item in trunk:
            b = bead.beadSequence(item)
            kmers = kmer.kmerCount(b, k)
            kmerNumber = len(kmers.kmers)
            kmerDist[kmerNumber] = kmerDist.get(kmerNumber, 0) + 1
        outputDist = sorted([(i, j) for i, j in kmerDist.items()], key=lambda x:x[0])
        print(outputDist[0:10])
        print(outputDist[-11:-1])

        f.write('kmerNumber\tCount\n')
        for line in outputDist:
            f.write('{0}\t{1}\n'.format(line[0], line[1]))