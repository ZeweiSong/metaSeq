#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu May 31 09:23:46 2018

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%%
from __future__ import print_function
from __future__ import division
from metaSeq import io as seqIO
import argparse
import json
import textwrap

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,\
                                 description=textwrap.dedent('''\
                                                             Remove replicated sequences from the assembled and unassembled stLFR data.
                                                             There is an option to output data as JSON format, which should be a faster
                                                             alternative.'''),\
                                     epilog=textwrap.dedent('''\
                                        ------------------------
                                        By Zewei Song
                                        Environmental Ecology Lab
                                        Institute of Metagenomics
                                        BGI-Research
                                        songzewei@genomics.cn
                                        songzewei@outlook.com
                                        ------------------------'''))
parser.add_argument('-single', help='The assembled file.')
parser.add_argument('-twin', help='The unassembled files, separated by comma: File1,File2')
parser.add_argument('-o', help='The output sequence file.')
parser.add_argument('-json', help='The JSON file to save the bead dictionary.')

args = parser.parse_args()

# Read in the assembled and unassembled sequences
# Assembled file
singleFile = args.single
seqSingle = seqIO.sequence(singleFile, fastx='q')

seqCount = 0

bead = {}
for record in seqSingle:
    if record[1].find('N') == -1: # Keep sequence without N
        barcode = record[0].split('/')[-1]
        try:
            bead[barcode]['s'][record[1]] += 1
        except KeyError:
            bead[barcode] = {'s':{}}
            bead[barcode]['s'][record[1]] = 1
        seqCount += 1
        if seqCount % 1000000 == 0:
            with open('bead_dereplicate.log', 'w') as f:
                f.write('Processed {0} sequences. Currently found {1} beads.'\
                        .format(seqCount, len(bead)))
    else:
        pass
    
# Unasembled forward and reverse reads
twinFile = args.twin.split(',')
seqTwin = seqIO.sequence_twin(twinFile[0], twinFile[1], fastx='q')

for r1, r2 in seqTwin:
    barcode = r1[0].split('/')[-1] # There is no check up for the correspondence of barcode in R1 and R2!
    seq = r1[1] + '_' + r2[1]
    if seq.find('N') == -1: # Keep sequence without N
        try:
            bead[barcode]['t'][seq] += 1
        except KeyError:
            try:
                bead[barcode]['t'] = {}
                bead[barcode]['t'][seq] = 1
            except KeyError:
                bead[barcode] = {'t':{}}
                bead[barcode]['t'][seq] = 1
        seqCount += 1
        if seqCount % 1000000 == 0:
            with open('bead_dereplicate.log', 'w') as f:
                f.write('Processed {0} sequences. Currently found {1} beads.'\
                        .format(seqCount, len(bead)))
    else:
        pass

with open('bead_dereplicate.log', 'w') as f:
    f.write('Processed {0} sequences. Found {1} beads.'.format(seqCount, len(bead)))

if args.json:
    with open(args.json, 'w') as f:
        json.dump(bead, f)
else:
    pass
#%%            
# Write into a FASTA file
outputFile = args.o
with open(outputFile, 'w') as f:
    for key, value in bead.items():
        count = 0
        try:
            for seq in value['s'].keys():
                label = key + '-' + str(count)
                f.write('>%s\n' % label)
                f.write('%s\n' % seq)
                count += 1
        except KeyError:
            pass
        
        try:
            for seqs in value['t'].keys():
                seqs = seqs.split('_')
                f.write('>{0}\n'.format(key + '-0-' + str(count)))
                f.write('{0}\n'.format(seqs[0]))
                f.write('>{0}\n'.format(key + '-1-' + str(count)))
                f.write('{0}\n'.format(seqs[1]))
                count += 1
        except KeyError:
            pass