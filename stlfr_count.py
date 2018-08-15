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
with open(inputFile, 'r') as f:
    for line in f:
        bead = json.loads(line.strip('\n'))