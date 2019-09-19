#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 10:36:47 2018

Coders who love to comment their code are unlikely to have bad luck.

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%%
from __future__ import print_function
from __future__ import division

''' Calculate the kemr number distribution
kmerDist = {}

for i in range(40):
    print(i)
    with open('kmerDist/dist_{0}.tsv'.format(i), 'r') as f:
        f.readline() # skip the header
        for line in f:
            line = line.strip('\n').split('\t')
            kmerDist[int(line[0])] = kmerDist.get(int(line[0]), 0) + int(int(line[1]))

kmers = sorted([(i, j) for i, j in kmerDist.items()])
print(kmers[-1])

with open('kmerDist/dist.tsv', 'w') as f:
    f.write('kmerNumber\tCount\n')
    for item in kmers:
        f.write('{0}\t{1}\n'.format(item[0], item[1]))
'''

''' Calculate the jacarrd distribution
bin = {} # Create ten bins for (0, 1]
for i in range(99):
    bin[(i/1000, (i+1)/1000)] = 0
bin[0] = 0
bin[0.01] = 0


def findBin(number):
    if number == 0:
        return 0
    elif number >= 0.01:
        return 0.01
    else:
        return (number//0.001/1000, (number//0.001+1)/1000)

#print(bin)
#print(findBin(0.0023))
#%%
content = []
with open('kmer.sample.pd.tsv', 'r') as f:
    print('Header is {0}'.format(f.readline()))
    for line in f:
        line = line.strip('\n').split('\t')
        line[0] = line[0][0:-1]
        line[2] = float(line[2])
        content.append(line)
#        binNumber = findBin(line[2])
#        bin[binNumber] += 1
print(len(content))
#%%
with open('kmer.jcd.0.02.tsv', 'w') as f:
    f.write('Target\tSource\tJacarrd\n')
    for item in content:
        if item[2] > 0.02:
            f.write('{0}\t{1}\t{2}\n'.format(item[0], item[1], item[2]))
'''

#%%
''' Extract bead sequences by module number '''
from metaSeq import io as seqIO
from metaSeq import bead

module = {}
with open('kmer.jcd.0.02.module.txt', 'r') as f:
    f.readline()
    for line in f:
        line = line.strip('\n').split('\t')
        module[line[0]] = line[1]
print(len(module))

cluster = {}

for item in list(set(module.values())):
    cluster[item] = []
print(len(cluster))        
beads = seqIO.beadJson('CL100077200_L01.json')
for item in beads:
    b = bead.beadSequence(item)
    classNumber = module.get(b.barcode, False)
    if classNumber:
        cluster[classNumber] += b.fastaSequences()
print(len(cluster))
for key, value in cluster.items():
    seqIO.write_seqs(value, 'cluster/{0}.fa'.format(key), fastx='a', mode='w')