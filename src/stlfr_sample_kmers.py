#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 10:36:47 2018

Coders who love to comment their code are unlikely to have bad luck.

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%%
from __future__ import print_function
from __future__ import division
from metaSeq import io as seqIO
import json
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='JSON kmer file.')
parser.add_argument('-o', '--output', help='Output file.')
parser.add_argument('-r', '--range', help='A range for sampling kmer by kmer number, use "[]" to include margin, and "()" to exclude margin.')
args = parser.parse_args()
args = parser.parse_args(['-i', 'kmer.json', '-o', 'kmer.sample.json', '-r', '[2300, 10000]'])
inputFile = args.input
outputFile = args.output
rang = args.range
left = int(rang[1:-1].split(',')[0])
right = int(rang[1:-1].split(',')[1])
if rang[0] == '[':
    left -= 1
else:
    pass
if rang[-1] == ']':
    right += 1
else:
    pass

count = 0
countPass = 0
beads = seqIO.beadJson(inputFile)
with open(outputFile, 'w') as f:
    for item in beads:
        count += 1
        kmerNumber = len(list(item.values())[0])
        if kmerNumber > left and kmerNumber < right:
            countPass += 1
            f.write('{0}\n'.format(json.dumps(item)))
        print(count, countPass)