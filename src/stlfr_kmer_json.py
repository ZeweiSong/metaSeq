#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 30 07:35:33 2019

@author: Zewei Song
@email: songzewei@genomics.cn
"""
from __future__ import print_function
from __future__ import division
import argparse
import json
from metaSeq import kmer
from metaSeq import bead

parser = argparse.ArgumentParser()
parser.add_argument('-i', help='Input JSON-Bead file.')
parser.add_argument('-o', help='Output JSON-Bead-Kmer file.')
parser.add_argument('-k', default=21, type=int, help='Kmer size')
parser.add_argument('-report', default='kmerReport.tsv', help='Report file on Kmers.')
args = parser.parse_args()

inputFile = args.i
outputFile = args.o
kmerSize = args.k
reportFile = args.report

kmerPool = []
report = {}
kmerFrag = []
with open(inputFile, 'r') as f:
    for line in f:
        b = bead.beadSequence(json.loads(line))
        kmers = kmer.kmerCount(b, kmerSize)
        barcode = b.barcode
        kmerPool.append({barcode:kmers.kmers})
        kmerNumber = len(kmers.kmers)
        fragNumber = len(b.fragments)
        kmerFrag.append((kmerNumber, fragNumber))
        report[kmerNumber] = report.get(kmerNumber, 1) + 1
with open(outputFile, 'w') as f:
    for line in kmerPool:
        f.write('{0}\n'.format(json.dumps(line)))

report = sorted([x for x in report.items()], key=lambda i:i[0])

with open(reportFile, 'w') as f:
    f.write('KmerNumber\tCount\n')
    for line in report:
        f.write('{0}\t{1}\n'.format(line[0], line[1]))

with open('kmerFrag.tsv', 'w') as f:
    f.write('KmerNumber\tFragNumber\n')
    for line in kmerFrag:
        f.write('{0}\t{1}\n'.format(line[0], line[1]))