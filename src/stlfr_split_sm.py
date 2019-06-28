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
from metaSeq import io as seqIO
from metaSeq import barcode as seqBar
import textwrap
import argparse
import os
import sys
import time
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
parser.add_argument('-bl', default='42', help='Length of the barcode tail on R2.')
parser.add_argument('-o', '--output', default='stlfrSplit', help='Name of the output folder.')
#parser.add_argument('-o', help='Output file, suffix will be added for four files.')
args = parser.parse_args()
#args = parser.parse_args(['-r1', 'r1.fq.gz', '-r2', 'r2.fq.gz', '-bl', '54', '-o', 'stlfrSplit'])
r1File = args.r1
r2File = args.r2
bl = int(args.bl)
outputFolder = args.output
if os.path.isdir(outputFolder):
    print('The folder "{0}" exist. Please remove the folder by hand.'.format(outputFolder))
    sys.exit()
else:
    os.mkdir(outputFolder)
t1 = time.time()
#outputFile = args.o

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
def barcode_set(seq, bl=bl, offset=0):
    if bl == 54:
        return [seq[100+offset:110+offset], seq[116+offset:126+offset], seq[144+offset:154+offset]]
    elif bl == 42:
        return [seq[100+offset:110+offset], seq[116+offset:126+offset], seq[132+offset:142+offset]]


# Return the number tuple of current read
def number_tuple(barcodes, NoSnpDict, OneSnpDict, trunkSize = 256): # barcodes is a list of three 10 bp string
    number = []
    # The majority of barcodes should be in the NoSnoDict
    for item in barcodes:
        try:
            number.append(NoSnpDict[item])
        except KeyError:
            try:
                number.append(OneSnpDict[item])
            except KeyError:
                return None
    return tuple(number)
    # return tuple((i//(trunkSize+1)+1 for i in number))


# Convert barcode number to trunk ordination
def number2ord(number_tuple):
    return tuple(i//257 + 1 for i in number_tuple)

#%% Read in the barcode list
# Forward and Reverse barcodes are saved in two Dictionaries.
print('Reading in the barcode list ...')
# Get the dict for no snp barcode
barcodeDictForward = seqBar.stlfrBarcode()
barcodeDictReverse = {}
NoSnpDict = {}

for key, value in barcodeDictForward.items():
    barcodeDictReverse[key] = rc(value[0])
for key, value in barcodeDictForward.items():
    for barcode in value:
        NoSnpDict[barcode] = int(key)
for key, value in barcodeDictReverse.items():
    for barcode in value:
        NoSnpDict[barcode] = int(key)

# Get the dict for 1 snp barcode
barcodeSnpDictForward = {}
barcodeSnpDictReverse = {}
OneSnpDict = {}

# Add all possible 1 SNP mutations for all barcode
barcodeSnpDictForward = {}
for key, value in barcodeDictForward.items():
    barcodeSnpDictForward[key] = snp_list(value[0])

# Create the RC Dict
barcodeSnpDictReverse = {}
for key, value in barcodeDictReverse.items():
    barcodeSnpDictReverse[key] = snp_list(value[0])

for key, value in barcodeSnpDictForward.items():
    for barcode in value:
        OneSnpDict[barcode] = int(key)
for key, value in barcodeSnpDictReverse.items():
    for barcode in value:
        OneSnpDict[barcode] = int(key)

print('No SNP dictionary has {0} keys.'.format(len(NoSnpDict)))
print('One SNP dictonary has {0} keys.'.format(len(OneSnpDict)))


#%% Read in R1 and R2 file
# The two number Dictionary contains all possible sequences with 1 base Error
print('Start reading {0} and {1} . . .'.format(r1File, r2File))
print('Splited reads are writing under the folder {0} . . .'.format(outputFolder))

# Create the list of output files. There are 6x6x6 files, plus one for error.
outputFileDict = {(0,0,0): outputFolder + '/' + '0_0_0.fq'}
for x in range(1,7):
    for y in range(1,7):
        for z in range(1,7):
            outputFileDict[(x,y,z)] = '{3}/{0}_{1}_{2}.gp.fq'.format(x,y,z,outputFolder)
   
# Assign read into bins using the first barcode as key
with open('stlfr_split_sm.log', 'w') as f1:
    count = 0
    error_count = 0
    seqs = seqIO.sequence_twin(r1File, r2File)
    for r1, r2 in seqs:
        count += 1
        if count // 1000000 >= 1:
            f1.write('\tProcessed {8.2f} M reads.'.format(count // 1000000))
        bead = number_tuple(barcode_set(r2[1]), NoSnpDict, OneSnpDict)
        if bead:
            b = '_'.join([str(i) for i in bead])
            r1[0] = r1[0][:-2] + '/' + b + '/1'
            r2[0] = r2[0][:-2] + '/' + b + '/2'
            r2[1] = r2[1][:100]
            r2[3] = r2[3][:100]
            seqIO.write_seqs([r1,r2], outputFileDict[number2ord(bead)], fastx='q', mode='a', gz=False)
        else:
            error_count += 1
            r1[0] = r1[0][:-2] + '/' + '0_0_0' + '/1'
            r2[0] = r2[0][:-2] + '/' + '0_0_0' + '/1'
            r2[1] = r2[1][:100]
            r2[3] = r2[3][:100]
            seqIO.write_seqs([r1,r2], outputFileDict[(0,0,0)], fastx='q', mode='a', gz=False)
t2 = time.time()
with open('stlfr_split_sm.log', 'a') as f2:
    f2.write('{0} total seqs, {1:8.2f} M.\n'.format(count, count//1000000))
    f2.write('{0} seqs do not have eligible barcode.\n'.format(error_count))
    t2 = time.time()
    f2.write('Used {0:8.0f} seconds ({1:8.2f} hrs) for processing the reads.\n'.format((t2 - t1), (t2 - t1)/3600))