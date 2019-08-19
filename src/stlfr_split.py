#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 2017-12-11
Last update on 2018-05-30

This script will split the raw stLFR data into individual bead. Short read with
the same barcode will be group together.

The script takes four inputw:
    1) raw R1 read
    2) raw R2 read with 54 bp barcode string
    3) Barcode file contains the 1526 barcode sequences and corresponding number
    4) Base string for output

The script spits five outputw:
    1) R1 read orderred by barcode
    2) R2 read orderred by barcode (54 bp string trimmed)
    3) R1 read without eligible barcode
    4) R2 read without eligible barcode
    5) A log file on the progress of splitting per 1M read, and a final report

Each barcode can tolerant 1 mutation. The three barcodes are identify from the
54 bp string using the pattern 10 + 6 + 10 + 18 + 10.

For a general BGISEQ PE100 lane with 600M reads, the split usually takes 8 hrs.

The current output is in FASTA format. However, this is not an efficient format
for stLFR data. A better alternative is to save as a Python dictionary using
bead (barcode) as the key. This option will be added later or you can parse the
data into Python dictionary yourself.

You will need to request super memory at the super computer for this task. A low
memory version is possible, but may take much longer time. We don't have plan
for a low memory version in the short run unless there is a urgent need. Contact
us if you need one.

@author: Zewei Song
@email: songzewei@genomics.cn or songzewei@outlook.com
"""
#%%
from __future__ import print_function
from __future__ import division
from metaSeq import io as seqIO
import sys
import textwrap
import argparse
import time
import json
#import gzip
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
                                        Split raw stLFR data into beads. Open the script to see a detail document.
                                        ------------------------
                                        By Zewei Song
                                        Environmental Ecology Lab
                                        Institute of Metagenomics
                                        BGI-Research
                                        songzewei@genomics.cn
                                        songzewei@outlook.com
                                        ------------------------'''))
parser.add_argument('-r1', help='Read1 file of the stLFR library.')
parser.add_argument('-r2', help='Read2 file of the stLFR library.')
parser.add_argument('-b', default='barcode.list', help='Forward barcode list, this is a text file.')
parser.add_argument('-o', help='Output base for FASTQ or JSON.')
parser.add_argument('-fastq', action='store_true', help='Turn on output to FASTQ, suffix will be added for four files.')
parser.add_argument('-json', action='store_true', help='Turn on output to JSON format, will write to one file with .json extension.')
parser.add_argument('-bl', default='42', help='Specify the length of barcode string, 42 or 54 bp.')
#parser.add_argument('-not_gz', action='store_false', help='Specify if the input is not a gz file, in most case you do not need it.')
args = parser.parse_args()
r1File = args.r1
r2File = args.r2
barcodeFile = args.b
base = args.o
bl = args.bl
not_gz = args.not_gz
t1 = time.time()

if not args.fastq and not args.json:
    print('Please specify at least one output format.\n')
    sys.exit()
#%% Functions
# Return the Reverse Compliment of a sequence
def rc(seq):
    replace_dict = {'A':'T', 'T':'A','C':'G','G':'C'}
    reverse = seq[::-1].upper()
    rc = ''
    for base in reverse:
        rc += replace_dict[base]
    return rc


# Generate all possible 1 SNP mutation of a given sequence
    # Original sequence is NOT included in the output
    # RC is not considered
    # There is l x 3 + 1 sequences including the original one
def snp_list(seq):
    snp_dict = {'A':['T','C','G'], 'T':['A','C','G'],\
                'C':['A','T','G'], 'G':['A','T','C']}
    snp_seqs = []
    for index, base in enumerate(seq): # mutate per base
        for mutate in snp_dict[base]:
            seqMutate = seq
            seqMutate = seqMutate[:index] + mutate + seqMutate[index+1:]
            snp_seqs.append(seqMutate)
    return snp_seqs


# Detect the three 10 bp barcode in the last 54 bp of Read 2
    # Currently using: 10 + 6 + 10 + 18 + 10 (three sets of 10 base barcode)
    # Lates barcode is 42bp 10 + 6 + 10 + 6 + 10 (Cuurently used by CNGB)
    # BAsed on our latest test on 2018-6-15, they may accidentally sequence R2
    # as 154 bp while using the 42 bp system. If this is the case, use the first
    # 142 bp.
def barcode_set(seq, bl, offset):
    if bl == '54':
        return [seq[100+offset:110+offset], seq[116+offset:126+offset], seq[144+offset:154+offset]]
    elif bl == '42':
        return [seq[100+offset:110+offset], seq[116+offset:126+offset], seq[132+offset:142+offset]]


# Return the number barcode set if exist
def number_set(barcodes, forwardDict, forwardSnpDict, reverseDict):
    number = []
    for item in barcodes:
        try:
            number.append(forwardDict[item])
        except:
            try:
                number.append(forwardSnpDict[item])
            except KeyError:
                #print('forward: {0}'.format(item))
                break
    if len(number) == 3: # barcode set (all three barcodes) found in the forward direction)
        return '_'.join(number)
    else: # Does not found barcode in the forward direction
        number = []
        for item in barcodes:
            try:
                number.append(reverseDict[rc(item)])
            except KeyError: # At least one barcode not found in reverse direction either
                #print('reverse: {0}'.format(item))
                return None
    return '_'.join(number)

''' Dyfunced now, RIP.
# Iterator for two files
# It only work for files with ABSOLUTELY corresponding record.
class sequence_twin(object):
    def __init__(self, file_r1, file_r2, fastx='a', gz=False):
        self.fastx = fastx
        self.gzip = gz
        if self.gzip:
            import gzip
            self.r1 = gzip.open(file_r1, 'rt')
            self.r2 = gzip.open(file_r2, 'rt')
        else:
            self.r1 = open(file_r1, 'r')
            self.r2 = open(file_r2, 'r')
        if fastx == 'a': self.n = 2
        elif fastx == 'q': self.n = 4
        else:
            print('Please specify the right format, "a" for FASTA and "q" for FASTQ.')
            self.n = 1

    def __iter__(self):
        return self

    def __next__(self):
        record = [[],[]]
        for i in range(self.n):
            line_r1 = self.r1.readline().strip('\n')
            line_r2 = self.r2.readline().strip('\n')
            if line_r1:
                record[0].append(line_r1)
                record[1].append(line_r2)
            else:
                raise StopIteration
        record[0][0] = record[0][0][1:]
        record[1][0] = record[1][0][1:]
        return record[0], record[1]
'''
#%% Read in the barcode list
# Forward and Reverse barcodes are saved in two Dictionaries.
print('Reading in the barcode list from {0} ...'.format(barcodeFile))
# Save to a Dict use barcode number as key and barcode sequence as value
barcodeDictForward = {}
with open(barcodeFile, 'r') as f:
    for line in f:
        line = line.strip('\n').split('\t')
        barcodeDictForward[line[1]] = [line[0]]

# Add all possible 1 SNP mutations for all barcode
barcodeSnpDictForward = {}
for key, value in barcodeDictForward.items():
    barcodeSnpDictForward[key] = snp_list(value[0])

# Create the RC Dict
barcodeDictReverse = {}
for key, value in barcodeDictForward.items():
    rc_seq = rc(value[0]) # The first sequence is the original Forward barcode
    barcodeDictReverse[key] = [rc_seq] + snp_list(rc_seq)

# Convert number:barcode to barcode:number
numberDictForward = {}
numberSnpDictForward = {}
numberDictReverse = {}

for key, value in barcodeDictForward.items():
    for barcode in value:
        numberDictForward[barcode] = key
for key, value in numberSnpDictForward.items():
    for barcode in value:
        numberSnpDictForward[barcode] = key
for key, value in barcodeDictReverse.items():
    for barcode in value:
        numberDictReverse[barcode] = key
print('{0} unique barcode possibilities in Forward.'.format(len(numberDictForward)))
print('{0} unique barcode possibilities in Reverse.'.format(len(numberDictReverse)))
#with open('split.log.txt', 'w') as f:
#    f.write('Built barcode dictionary\n')


#%% Read in R1 and R2 file
# Assign barocde on the fly, and save in memory (This will take up a LARGE trunk of memory)
# The two number Dictionary contains all possible sequences with 1 base Error
beadDict = {} # This is the dictionary that organizes seqs by bead (barcode set)
beadError = {'0_0_0':[]}  # This is the dicionary that organizes seqs without barcode

logFile = base + '.log'

seqs = seqIO.sequence_twin(r1File, r2File)
count = 0
error_count = 0
offsets = [0,-1,1,-2,2]
for r1, r2 in seqs:
    count += 1
    if count % 1000000 == 0: # Report per 1 million reads
        with open(logFile, 'w') as f:
            f.write('Processed {0:8.2f} M reads\n'.format(count/1000000))
        #break
    for offset in offsets:
        bead = number_set(barcode_set(r2[1], bl,offset), numberDictForward, numberSnpDictForward, numberDictReverse)
        if bead:
            break

    if bead:
        r1[0] = r1[0][:-2] + '/' + bead + '/1'
        r2[0] = r2[0][:-2] + '/' + bead + '/2'
        r2[1] = r2[1][:100]
        r2[3] = r2[3][:100]
        try:
            beadDict[bead].append([r1, r2])
        except KeyError:
            beadDict[bead] = [[r1, r2]]
    else:
        r1[0] = r1[0][:-2] + '/' + '0_0_0' + '/1'
        r2[0] = r2[0][:-2] + '/' + '0_0_0' + '/1'
        r2[1] = r2[1][:100]
        r2[3] = r2[3][:100]
        beadError['0_0_0'].append([r1, r2])
        error_count += 1

with open(logFile, 'w') as f1:
    f1.write('Processed {0:8.2f} M reads\n'.format(count/1000000))
    f1.write('Processed finished.\n')
    f1.write('Found {0} beads.\n'.format(len(beadDict)))

    if args.fastq:
        r1OutputFile = base + '.r1_split.fq'
        r2OutputFile = base + '.r2_split.fq'
        r1Error = base + '.r1_error.fq'
        r2Error = base + '.r2_error.fq'

        print('Writing to R1 file ...')
        with open(r1OutputFile, 'w') as f:
            for key, value in beadDict.items():
                for record in value:
                    record[0][0] = '@' + record[0][0]
                    for line in record[0]:
                        f.write('{0}\n'.format(line))

        print('Writing to R2 file ...')
        with open(r2OutputFile, 'w') as f:
            for key, value in beadDict.items():
                for record in value:
                    record[1][0] = '@' + record[1][0]
                    for line in record[1]:
                        f.write('{0}\n'.format(line))

        print('Writing to R1 Error file ...')
        with open(r1Error, 'w') as f:
            for key, value in beadError.items():
                for record in value:
                    record[0][0] = '@' + record[0][0]
                    for line in record[0]:
                        f.write('{0}\n'.format(line))

        print('Writing to R2 Error file ...')
        with open(r2Error, 'w') as f:
            for key, value in beadError.items():
                for record in value:
                    record[1][0] = '@' + record[1][0]
                    for line in record[1]:
                        f.write('{0}\n'.format(line))
    if args.json:
        f1.write('Output as JSON format.\n' )
        jsonOutputFile = base + '.json'
        with open(jsonOutputFile, 'w') as f:
            for key, value in beadDict.items():
                currentBead = {key:[]}
                currentBead[key] = [seqIO.fastq2list(i[0], i[1]) for i in value]
                f.write('%s\n' % json.dumps(currentBead))

    f1.write('{0} total seqs.\n'.format(count))
    f1.write('{0} seqs do not have eligible barcode.\n'.format(error_count))
    t2 = time.time()
    f1.write('Used {0:8.0f} seconds ({1:8.2f} hrs) for processing the reads.\n'.format((t2 - t1), (t2 - t1)/3600))
