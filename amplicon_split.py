#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 11:42:28 2018

A place holder for split function.

@author: Zewei Song
@email: songzewei@genomics.cn
"""
from __future__ import print_function
from __future__ import division
import argparse
from metaSeq import io as seqIO

parser = argparse.ArgumentParser()
parser.add_argument('-r1', help='R1 reads')
parser.add_argument('-r2', help='R2 reads if needed')
parser.add_argument('-b', '--barcode', help='Barcode mapping file with sample info')
parser.add_argument('-o', '--output', help='Output folder for splitted files')
seqFormat = parser.add_mutually_exclusive_group()
seqFormat.add_argument('-gz', action='store_true')
seqFormat.add_argument('-fq', action='store_true')
args = parser.parse_args()

compress = True
if args.fq:
    compress = False

r1File = args.r1
if args.r2:
    r2File = args.r2
    seq = seqIO.sequence_twin(r1File, r2File, fastx='q', gz=compress)
else:
    seq = seqIO.sequence(r1File, fastx='q', gz=compress)

# Parse barcode list
barcode = {}
with open(args.barcode, 'r') as f:
    for line in f:
        line = line.strip('\n').split('\t')
        barcode[line[0]] = line[1]

# Get ready the output file name
# If R1 only, output as SampleName.r1.fq
# If R1 and R2, output as SampleName.r1.fq, SampleName.r2.fq
output = {}
for key, value in barcode.items():
    if args.r2:
        output[key] = [value + '.r1.fq', value + '.r2.fq']
    else:
        output[key] = [value + '.r1.fq']

# Get the barcode
# Barcode locates at the first 6 bp on both read
for record in seq:
    if args.r2:
        read1 = record[0]
        read2 = record[1]
        tag = [read1[1][:6], read2[1][:6]]
        if tag[0] == tag[1]:
            tag = tag[0]
            try:
                seqIO.write_seqs(read1, output[tag][0], fastx='q', mode='a')
                seqIO.write_seqs(read2, output[tag][1], fastx='q', mode='a')
            except KeyError:
                pass
    else:
        read1 = record
        tag = read1[1][:6]
        try:
            seqIO.write_seqs(read1, output[tag][0], fastx='q', mode='a')
        except KeyError:
            pass