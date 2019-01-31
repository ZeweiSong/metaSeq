#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 10:22:53 2019

Test on multiprocessing. If you see this sentence, I'm still trying. Stay tune ...

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%%
from __future__ import print_function
from __future__ import division
from multiprocessing import Pool
import json
from itertools import combinations
import argparse
from metaSeq import kmer
from metaSeq import bead

parser = argparse.ArgumentParser()
parser.add_argument('-i', help='Input JSON-Bead file.')
parser.add_argument('-t', '--threshold', nargs='+', default=[0.02,0.04], type = float, help='Threshold for distance.')
parser.add_argument('-rawout', help='Output raw edge file.')
parser.add_argument('-tout', help='Output thresholded edge file.')
parser.add_argument('-k', default=21, type=int, help='Kmer size default = 21')
parser.add_argument('-j', default = 10, type=int, help='Number of jobs to parallelize.')
args = parser.parse_args()

inputFile = args.i
outputRaw = args.rawout
outputThreshold = args.tout
threshold = args.threshold
kmerSize = args.k
job = args.j

# Read in JSON-Bead file
#Calculate kmer pools for all beads
kmerPool = []
beadCount = 0
with open(inputFile, 'r') as f:
    for line in f:
        b = bead.beadSequence(json.loads(line))
        kmerPool.append(kmer.kmerCount(b, kmerSize))
        beadCount += 1
print('Found {0} beads.'.format(beadCount))

'''
# Index pool for all index pairs
indexPool = []
for pair in combinations(list(range(len(kmerPool))),2):
    indexPool.append(pair)
print('We have to go over {0} pair for kmer distance.'.format(len(indexPool)))
'''

# Divide the kmer pair pool
pairPool = []
for pair in combinations(kmerPool,2):
    pairPool.append(pair)
#print(pairPool)
step = len(pairPool) // job
print('Step is {0}'.format(step))
start = 0
pairPoolPool = []
for i in range(job):
    if i+1 < job: # not the last job
        pairPoolPool.append(pairPool[start:step])
        start += step
        print('Start change to {0}'.format(start))
    else:
        pairPoolPool.append(pairPool[start:])
print('The kmer pool is divided into {0} sections.'.format(len(pairPoolPool)))

#%% Paralle each index pool
def kmerDistanceWorker(x):
    print('Worker started.')
    #print(x)
    d = []
    n = 0
    for kmerPair in x:
        D = kmer.kmerDistance((kmerPair[0].set, kmerPair[1].set)).mashDistance()
        d.append((kmerPair[0].barcode, kmerPair[1].barcode, D))
        n += 1
        if n %100000 == 0:
            print(n)
    return d
if __name__ == '__main__':
    p = Pool(job)
    result = p.map(kmerDistanceWorker, pairPoolPool)
    with open(outputRaw, 'w') as f:
        f.write('Source\tTarget\tDistance\n')
        for item in result:
            for line in item:
                f.write('{0}\t{1}\t{2}\n'.format(line[0], line[1], line[2]))