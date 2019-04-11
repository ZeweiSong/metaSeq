#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  6 10:47:32 2018

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%%
from __future__ import print_function
from __future__ import division
#from metaSeq import io as seqIO
import textwrap
import argparse
import time
import os
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
parser.add_argument('-o', help='Output file, suffix will be added for four files.')
parser.add_argument('-s', default=1000000, help='Number of record read in every time, default is 1 million.')
parser.add_argument('-not_gz', action='store_false', help='Specify if the input is gz file, in most case you do not need it.')
args = parser.parse_args()
r1File = args.r1
r2File = args.r2
barcodeFile = args.b
outputFile = args.o
not_gz = args.not_gz
t1 = time.time()

size = int(args.s)

#barcodeFile = 'barcode.list'
#r1File = 'r1.fq'
#r2File = 'r2.fq'
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
def barcode_set(seq):
    return [seq[100:110], seq[116:126], seq[144:154]]


# Return the number barcode set if exist
def number_set(barcodes, forwardDict, reverseDict):
    number = []
    for item in barcodes:
        try:
            number.append(forwardDict[item])
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


# Iterator for two files
# It only work for files with ABSOLUTELY corresponding record.
class sequence_twin_trunk(object):
    def __init__(self, file_r1, file_r2, fastx='a', gz=False, trunk_size=1000000):
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
        self.trunk_size = trunk_size

    def __iter__(self):
        return self
    
    def __next__(self):
        r1_trunk = []
        r2_trunk = []
        for record in range(self.trunk_size):
            r1 = []
            r2 = []
            for i in range(self.n):
                line_r1 = self.r1.readline().strip('\n')
                line_r2 = self.r2.readline().strip('\n')
                if line_r1:
                    r1.append(line_r1)
                    r2.append(line_r2)
                else:
                    if len(r1_trunk) > 0:
                        return r1_trunk, r2_trunk
                    else:
                        raise StopIteration
            r1[0] = r1[0][1:]
            r2[0] = r2[0][1:]
            r1_trunk.append(r1)
            r2_trunk.append(r2)
        return r1_trunk, r2_trunk

tempDict = {}
for i in range(1536):
    tempDict[i+1] = str(i+1) + '.temp'

#%% Save to a Dict use barcode number as key and barcode sequence as value
barcodeDictForward = {}
with open(barcodeFile, 'r') as f:
    for line in f:
        line = line.strip('\n').split('\t')
        barcodeDictForward[line[1]] = [line[0]]

# Add all possible 1 SNP mutations for all barcode
for key, value in barcodeDictForward.items():
    barcodeDictForward[key] += snp_list(value[0])

# Create the RC Dict
barcodeDictReverse = {}
for key, value in barcodeDictForward.items():
    rc_seq = rc(value[0]) # The first sequence is the original Forward barcode
    barcodeDictReverse[key] = [rc_seq] + snp_list(rc_seq)

# Convert number:barcode to barcode:number
numberDictForward = {}
numberDictReverse = {}

for key, value in barcodeDictForward.items():
    for barcode in value:
        numberDictForward[barcode] = key
for key, value in barcodeDictReverse.items():
    for barcode in value:
        numberDictReverse[barcode] = key


#%% Read in R1 and R2 file
# The two number Dictionary contains all possible sequences with 1 base Error

beadError = {'0_0_0':[]}  # This is the dicionary that organizes seqs without barcode
trunks = sequence_twin_trunk(r1File, r2File, fastx='q', trunk_size = size, gz=not_gz)

# Assign read into bins using the first barcode as key
for t1, t2 in trunks:
    beadDict = {} # This is the dictionary that organizes seqs by bead (barcode set)
    for i in range(1536):
        beadDict[i+1] = []
    for r1, r2 in zip(t1, t2):
        bead = number_set(barcode_set(r2[1]), numberDictForward, numberDictReverse)
        if bead:
            firstBarcode = int(bead.split('_')[0])
            beadDict[firstBarcode].append((bead, r1[1], r1[3], r2[1][:100], r2[3][:100]))
        else:
            beadError['0_0_0'].append((r1[1], r1[3], r2[1][:100], r2[3][:100]))

    # Write to temp file:
    for key, value in beadDict.items():
        if len(value) > 0:
            with open(tempDict[key], 'a') as f:
                for line in value:
                    f.write('{0}\n'.format('\t'.join(line)))

#%% Order each temp file by bead, and write to a new file
#outputFile = 'orderred.fq'
for key, value in tempDict.items():
    temp = []
    try:
        with open(value, 'r') as f:
            for line in f:
                temp.append(line.strip('\n').split('\t'))
        temp.sort(key=lambda x:x[0])
        with open(outputFile, 'a') as f:
            for line in temp:
                f.write('{0}\n'.format('\t'.join(line)))
    except FileNotFoundError:
        pass

#%% Remove temp files
for key, value in tempDict.items():
    if os.path.isfile(value):
        os.remove(value)