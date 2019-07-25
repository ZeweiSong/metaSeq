#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 10:36:47 2018

This is a demonstration on how we can use multiple alignment results.

Or to be precise, the alignment space with ALL eligible hits, not just the first hit.

We simulate genomes by randomly drawing from a component pool. Component is stable inside, 
and is totally different from other components.

Currently, there is no evolutionary relationship among genomes.

By given an array of abundance, we can simulate the total DNA pool, and randomly draw 
a shotgun pool.

--
Coders who love to comment their code are unlikely to have bad luck.

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%%
from __future__ import print_function
from __future__ import division
import numpy as np
import networkx as nx

# Function to create a genome from the parts pool with length l
def create(l, pool):
    return tuple(np.random.choice(pool, l, replace=False))

def intersection(lst1, lst2):
    temp = set(lst2)
    lst3 = [value for value in lst1 if value in temp]
    return lst3

# Set the size of componet pool
size = 100
parts = [i for i in range(size)] # A pool for genomic basic units

# Species pool, genome is genome, abundance is abundance
species = {}
# s = (('a',3,20),('b',4,10),('c',2,40),('d',3,90),('e',5,70)) # Abundance array for (name, genome length, abundance)
s = (('a',80,1000),('b',80,1000)) # A two species sample
for item in s:
    species[item[0]] = {}
    species[item[0]]['genome'] = create(item[1], parts)
    species[item[0]]['abundance'] = item[2]

print('A is {0}, B is {1}.'.format(len(species['a']['genome']), len(species['b']['genome'])))
shared = intersection(species['a']['genome'], species['b']['genome'])
print('A and B share {0}.'.format(len(shared)))

# Get the overlap between A and B

#%%
# Create an index from component to genomes
index = {}
for key, value in species.items():
    for item in value['genome']:
        index[item] = index.get(item, []) + [key]
print('The index length is {0}.'.format(len(index)))

# Get The total DNA pool. this is all the components we can find in this sample, ALL OF'EM
dna = []
for key, value in species.items():
    dna += value['genome'] * value['abundance']
print('The DNA pool contians a total of {0} components.'.format(len(dna)))


#%% Simulated sequencing
# Draw a shallow shot gun pool
depth = len(dna)//100 # sequencing depth, always lower than the size of DNA pool
shotgun = np.random.choice(dna, depth, replace=False)
print('The shotgun seq results in {0} readed components.'.format(len(shotgun)))
print('Depth is {0}.'.format(len(shotgun)/len(dna)))
#%% Simulated alignments
# Create the alignment space (all eligible alignments)
# Use index number as the name of each query
aln = []
for i, item in enumerate(shotgun):
    for ref in index[item]:
        aln.append((i, ref))
print('The shotgun data results in {0} eligible alignments to A and B.'.format(len(aln)))

#%% Dereplicate alignments
alnDerep = {}
alnDict = {}

for item in aln:
    alnDict[item[0]] = alnDict.get(item[0], []) + [item[1]]
for key in alnDict.keys():
    alnDict[key] = tuple(sorted(alnDict[key]))
for item in alnDict.values():
    alnDerep[item] = alnDerep.get(item, 0) + 1
print(alnDerep)

#%%
# Features of the alignment space
# Get the alignment space distribution
distCount = {}
for item in aln:
    distCount[item[0]] = distCount.get(item[0], 0) + 1
#print(distCount)

dist = {}
for key, value in distCount.items():
    dist[value] = dist.get(value, 0) + 1
print(dist)
#%%
# Create the QR-graph
# Assign attribute to query and ref nodes
G = nx.Graph()
for item in aln:
    G.add_node(item[0], query=True, wt=1)
    G.add_node(item[1], ref=True, wt=0)
    G.add_edge(item[0], item[1])

#print(G.nodes())
#print(G.edges())

# Use attribute to get the query and ref nodes
queryNodes = nx.get_node_attributes(G, "query").keys()
refNodes = nx.get_node_attributes(G, "ref").keys()
print('Found {0} queries.'.format(len(queryNodes)))
print('Found {0} refs.'.format(len(refNodes)))

# Find all subgraph of the graph. These subgraphs are not connected. Thus are 
# the first step solution for the alignment space
sg = [] # Create a list of subgraphs

while G.nodes():
    nodes = list(G.nodes().keys())
    nodes_to_remove = nx.shortest_path(G, nodes[0]).keys()
    sg.append(G.subgraph(nodes_to_remove).copy())
    G.remove_nodes_from(nodes_to_remove)

print('Graph divided into {0} subgraphs.'.format(len(sg)))
for item in sg:
    refNodes = nx.get_node_attributes(item, 'ref').keys()
    print(refNodes)

#%%
# New WTA use unique alignment as reference weight
# STEP 1, assign weight to all ref nodes.
queryNodes = nx.get_node_attributes(sg[0], 'query').keys() # Get all query node
for node in queryNodes:
    if sg[0].degree(node) == 1:
        refNode = list(sg[0].neighbors(node))[0]
        sg[0].node[refNode]['wt'] += 1
print(sg[0].node['a'])
print(sg[0].node['b'])

# STEP 2, allocate all tied query nodes among refs using their weights.
queryNodes = nx.get_node_attributes(sg[0], 'query').keys() # Get all query node
for node in queryNodes:
    d = {}
    if sg[0].degree(node) > 1:
        for item in sg[0].neighbors(node):
            d[item] = sg[0].node[item]['wt']
    refSum = sum([i for i in d.values()])
    for item in d.keys():
        sg[0].node[item]['wt'] += sg[0].node[node]['wt'] * (d[item] / refSum)
# STEP 3, return the abudnance of this subgraph
print(sg[0].node['a'])
print(sg[0].node['b'])