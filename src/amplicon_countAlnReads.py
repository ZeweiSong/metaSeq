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
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='Alignment to count.')
args = parser.parse_args()
inputFile = args.input

reads = {}
with open(inputFile, 'r') as f:
    for line in f:
        line = line.strip('\n').split('\t')
        currentRead = line[0]
        try:
            reads[currentRead] += 1
        except KeyError:
            reads[currentRead] = 1

# print('There are {0} reads in the alignment {1}.'.format(len(reads), inputFile))
print(len(reads))