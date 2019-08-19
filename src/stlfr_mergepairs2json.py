#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  8 12:44:11 2018

@author: Zewei Song
@email: songzewei@genomics.cn
"""
from __future__ import print_function
from __future__ import division
import argparse
from metaSeq import io as seqIO
import json
import textwrap

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
                                        Convert the pair end merged FASTQ files into JSON format.
                                        The JSON file can be read in using the Class bead_json().
                                        By Zewei Song
                                        Environmental Ecology Lab
                                        Institute of Metagenomics
                                        BGI-Research
                                        songzewei@genomics.cn
                                        songzewei@outlook.com
                                        ------------------------'''))
parser.add_argument('-a', help='FASTQ file of assembled reads.')
parser.add_argument('-f', help='Unassembled forward reads.')
parser.add_argument('-r', help='Unassembled reverse reads.')
parser.add_argument('-o', help='Output bead file in JSON format.')
#parser.add_argument('-fasta', action='store_true', default=False)
args = parser.parse_args()
assemFile = args.a
fwdFile = args.f
revFile = args.r
outputFile = args.o
logFile = 'stlfr_mergepairs2json.log'
#if args.fasta:
#    fastx = 'a'
#else:
#    fastx = 'q'
with open(logFile, 'w') as f:
    f.write('Reading in the files.\n')

bead = {} # Dictionary for storing bead
bead = seqIO.mergepairs2bead(assemFile, fwdFile, revFile)

with open(logFile, 'a') as f:
    f.write('Writing to JSON file.\n')
with open(outputFile, 'w') as f:
    for key, value in bead.items():
        f.write('%s\n' % json.dumps({key:value}))