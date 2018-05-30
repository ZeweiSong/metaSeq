#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 16:59:32 2017

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%%
from __future__ import print_function
from __future__ import division
from metaSeq import io
import argparse
import time
parser = argparse.ArgumentParser()
parser.add_argument('-r1', help='Read1 file')
parser.add_argument('-r2', help='Read2 file')
parser.add_argument('-b', help='Forward barcode list')
parser.add_argument('-o', help='Output base')
parser.add_argument('-not_gz', action='store_false', help='Specify if the input is gz file')
args = parser.parse_args()
r1File = args.r1
r2File = args.r2
barcodeFile = args.b
base = args.o
not_gz = args.not_gz

t1 = time.time()
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

seqs = io.sequence_twin(r1File, r2File, fastx='q', gz=not_gz)
count = 0
error_count = 0
for r1, r2 in seqs:
    count += 1
    if count % 1000000 == 0: # Report per 1 million reads
        with open(logFile, 'w') as f:
            f.write('Processed {0:8.2f} M reads\n'.format(count/1000000))
    bead = number_set(barcode_set(r2[1]), numberDictForward, numberDictReverse)
    if bead:
        r1[0] = r1[0] + '/' + bead
        r2[0] = r2[0] + '/' + bead
        r2[1] = r2[1][:100]
        r2[3] = r2[3][:100]
        try:
            beadDict[bead].append([r1, r2])
        except KeyError:
            beadDict[bead] = [[r1, r2]]
    else:
        r1[0] = r1[0] + '/' + '0_0_0'
        r2[0] = r2[0] + '/' + '0_0_0'
        r2[1] = r2[1][:100]
        r2[3] = r2[3][:100]
        beadError['0_0_0'].append([r1, r2])
        error_count += 1

with open(logFile, 'w') as f1:
    f1.write('Processed {0:8.2f} M reads\n'.format(count/1000000))
    f1.write('Processed finished.\n')
    f1.write('Found {0} beads.\n'.format(len(beadDict)))
    
    r1OutputFile = base + '.r1_split.fq'
    r2OutputFile = base + '.r2.split.fq'
    r1Error = base + '.r1_error.fq'
    r2Error = base + '.r2_error.fq'
    
    print('Writing to R1 file ...')
    with open(r1OutputFile, 'w') as f:
        for key, value in beadDict.items():
            for record in value:
                for line in record[0]:
                    f.write('{0}\n'.format(line))
    
    print('Writing to R2 file ...')
    with open(r2OutputFile, 'w') as f:
        for key, value in beadDict.items():
            for record in value:
                for line in record[1]:
                    f.write('{0}\n'.format(line))
    
    print('Writing to R1 Error file ...')
    with open(r1Error, 'w') as f:
        for key, value in beadError.items():
            for record in value:
                for line in record[0]:
                    f.write('{0}\n'.format(line))
    
    print('Writing to R2 Error file ...')
    with open(r2Error, 'w') as f:
        for key, value in beadError.items():
            for record in value:
                for line in record[1]:
                    f.write('{0}\n'.format(line))
    
    f1.write('{0} seqs have eligible barcode.\n'.format(count))
    f1.write('{0} seqs do not have eligible barcode.\n'.format(error_count))
    t2 = time.time()
    f1.write('Used {0:8.0f} seconds ({1:8.2f} hrs) for processing the reads..\n'.format((t2 - t1), (t2 - t1)/3600))