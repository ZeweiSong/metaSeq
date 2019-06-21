#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 10:36:47 2018

Count the number of reads that have at least one alignment found in the reference database.

Coders who love to comment their code are unlikely to have bad luck.

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%
from __future__ import print_function
from __future__ import division
import argparse
from metaSeq import io as seqIO
from metaSeq import amplicon
parser = argparse.ArgumentParser()
parser.add_argument('-d', '--direction', help='File after rule directin.')
parser.add_argument('-q', '--quality', help='File after rule quality.')
parser.add_argument('-a', '--align', help='File after rule align and rule keep_paired.')
args = parser.parse_args()
directionFile = args.direction
qualityFile = args.quality
alignFile = args.align

count = {'d':0, 'q':0, 'a': amplicon.countAlignment(alignFile)}
for item in seqIO.sequence(directionFile):
    count['d'] += 1
for item in seqIO.sequence(qualityFile):
    count['q'] += 1

print('Direction\t{0}\tQuality\t{1}\tAlignment\t{2}\n'.format(count['d'], count['q'], count['a']))