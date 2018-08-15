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
from metaSeq import bead
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', help='Input file to sample')
parser.add_argument('-o', help='Output file')
parser.add_argument('-n', help='Number of first n bead to output')
args = parser.parse_args()

inputFile = args.i
outputFile = args.o
n = args.n

beadParser = seqIO.stlfr_bead(inputFile)
for item in beadParser:
    seqIO.write_seqs(item[1], outputFile, fastx='q', mode='a')