#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 16:59:32 2017

@author: songz
"""
from __future__ import print_function
from __future__ import division
from metaSeq import io
from metaSeq import qc
import argparse
import time

parser = argparse.ArgumentParser()
parser.add_argument('-i', help='Input FASTQ file.')
parser.add_argument('-o', help='Output FASTQ file.')
parser.add_argument('-r', help='Expected Error rate, EE/bp')
parser.add_argument('-l', help='Minimum length to keep, >=')
parser.add_argument('-t', default=100000, help='Trunk size. By default is 100,000')
args = parser.parse_args()
input_file = args.i
output_file = args.o
rate = float(args.r)
ml = int(args.l)
trunk = int(args.t)
t1 = time.clock()
print('Reading {0} by trunk ({1}) ...'.format(input_file, trunk))
#%% Truncate filter at maxE rate 0.01
input_seq = io.sequence_trunk(input_file, fastx='q', trunk_size=trunk)
p_dict = qc.qual_score()

i = 0
j = 0
c = 0
for trunk in input_seq:
    c += 1
    content = []
    for record in trunk:
        i += 1
        filtered = qc.trunc_ee_rate(record, p_dict, rate=rate)
        if len(filtered[1]) >= ml:
            content.append(filtered)
            j += 1
    count = io.write_seqs(content, output_file, fastx='q', mode='a')
print('Expected error rate = {0}, Minimum length = {1}'.format(rate, ml))
print('Filtered {0}, kept {1}'.format(i, j))
t2 = time.clock()
print('Use {0} s.'.format(t2-t1))