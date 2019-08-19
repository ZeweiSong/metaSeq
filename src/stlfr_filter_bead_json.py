#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 30 17:33:41 2019

@author: Zewei Song
@email: songzewei@genomics.cn
"""
from __future__ import print_function
from __future__ import division
import argparse
import json
from metaSeq import bead

parser = argparse.ArgumentParser()
parser.add_argument('-i', help='Input JSON-Bead file.')
parser.add_argument('-o', help='Output JSON-Bead file.')
parser.add_argument('-fr', nargs='+', type=int, help='Fragment range, enter two int number MIN MAX (included)')
args = parser.parse_args()

inputFile = args.i
outputFile = args.o
minFrag = args.fr[0]
maxFrag = args.fr[1]

nKeep = 0
nDiscard = 0
with open(inputFile, 'r') as f:
    with open(outputFile, 'w') as f2:
        for line in f:
            b = bead.beadSequence(json.loads(line))
            fragNumber = len(b.fragments)
            if minFrag <= fragNumber <= maxFrag:
                f2.write('{0}\n'.format(json.dumps(b.json)))
                nKeep += 1
            else:
                nDiscard += 1

print('Within range {0} and {1}, {2} kept, {3} discarded.'.format(minFrag, maxFrag, nKeep, nDiscard))