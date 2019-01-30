#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 09:05:01 2019

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%%
from __future__ import print_function
from __future__ import division
import json
import networkx as nx
from itertools import combinations
from metaSeq import bead
from metaSeq import kmer
'''
inputJson = 'mock.filter.json'
kmerSize = 21

### Calculate the kmers of each bead, and write to a new json file
### Get the non-redundant kmer pool for the dataset
kmerJson = {}
n = 0
with open(inputJson, 'r') as f:
    for line in f:
        currentBead = bead.beadSequence(json.loads(line)) # Create a bead class from the JSON line
        barcode = currentBead.barcode
        print(barcode)
        n += 1
        print(n)
        kmers = kmer.kmerCount(currentBead, kmerSize)
        for item in kmers.kmers:
            #kmerJson[item] = kmerJson.get(item, [barcode]) + [barcode] # Causing additional barcode?
            try:
                kmerJson[item].append(barcode)
            except KeyError:
                kmerJson[item] = [barcode]

print(len(kmerJson))

with open('kmerJson2.json', 'w') as f:
    for key, value in kmerJson.items():
        f.write('%s\n' % json.dumps({key:value}))
print('Wrote to JSON')
# kmerJson is the dereplicated kmer table for the dataset
# The vale of the kmer is all associated bead
'''
### Get the bead graph using shared kmer as count
G = nx.Graph()
n = 0
with open('kmerJson2.json', 'r') as f:
    for line in f:
        k = json.loads(line)
        key = list(k.keys())[0]
        n += 1
        for pair in combinations(k[key], 2):
            try:
                G[pair[0]][pair[1]]['count'] += 1
            except KeyError:
                G.add_edge(pair[0], pair[1], count = 1)
        print(n, G.number_of_nodes(), G.number_of_edges())

### Add the number of kmer to each node


