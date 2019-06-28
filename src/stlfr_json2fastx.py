#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 11:25:52 2018

@author: Zewei Song
@email: songzewei@genomics.cn
"""
from __future__ import print_function
from __future__ import division
import argparse
from metaSeq import io as seqIO
from metaSeq import bead
import json

parser = argparse.ArgumentParser()
parser.add_argument('-i', help='Input beadJson format file.')
parser.add_argument('-o', help='The output FASTQ file base. Barcode will be put at -2, while /1 and /2 at -1.')
parser.add_argument('-bp', default=-1, type=int, help='The position of barcode in beadJson. usually -1 or -2')
parser.add_argument('-fastq', action='store_true', default=False, help='Enable to output FASTQ format.')
args = parser.parse_args()

inputFile = args.i
outputFileBase = args.o
barcodePosition = args.bp
fastq = args.fastq

if fastq:
    with open(inputFile, 'r') as f:
        for line in f:
            currentBead = bead.beadSequence(json.loads(line.strip('\n')))
            outputR1 = []
            outputR2 = []
            for record in currentBead.fragments:
                r1 = record[0:4]
                r2 = record[4:8]
                label = r1[0].split('/')
                r1[0] = label[0] + '/' + label[2] + '/' + label[1]
                label = r2[0].split('/')
                r2[0] = label[0] + '/' + label[2] + '/' + label[1]
                outputR1.append(r1)
                outputR2.append(r2)
            count = seqIO.write_seqs(outputR1, outputFileBase + '.r1.fq', fastx='q', mode='a')
            count = seqIO.write_seqs(outputR2, outputFileBase + '.r2.fq', fastx='q', mode='a')

else:
    with open(inputFile, 'r') as f:
        for line in f:
            currentBead = bead.beadSequence(json.loads(line.strip('\n')))
            outputFasta = currentBead.fastaSequences()
            count = seqIO.write_seqs(outputFasta, outputFileBase + '.fasta', fastx='a', mode='a')
    print('Wrote {0} sequences to {1}.'.format(count, outputFileBase + '.fasta'))