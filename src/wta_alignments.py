#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 10:36:47 2018

Analysis a full alignment space (alignments with all eligible hits)

Coders who love to comment their code are unlikely to have bad luck.

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%%
from __future__ import print_function
from __future__ import division
import networkx as nx

# load in alignment
alnFile = 's551.silva.b6'
alnOutput = 's551.silva.tsv'
aln = []
with open(alnFile, 'r') as f:
    for line in f:
        line = line.strip('\n').split('\t')
        aln.append((line[0], line[1]))
print('found {0} alignments.'.format(len(aln)))

'''
with open(alnOutput, 'w') as f:
    f.write('Target\tSource\n')
    for line in aln:
        f.write('{0}\t{1}\n'.format(line[0], line[1]))
'''
#%% Dereplicate the alignment
alnDerep = {}
alnDict = {}

for item in aln:
    alnDict[item[0]] = alnDict.get(item[0], []) + [item[1]] # item[0] is the query, item[1] is the reference
for key in alnDict.keys():
    alnDict[key] = tuple(sorted(alnDict[key])) # Sort the reference array
for item in alnDict.values(): # Count duplicated reference array
    alnDerep[item] = alnDerep.get(item, 0) + 1

# Filter alignment to keep Query <= k hits
k = 10
alnDerepFilter = {}
for key, value in alnDerep.items():
    if len(key) <= 10:
        alnDerepFilter[key] = value
alnSorted = sorted([(key, value) for key, value in alnDerepFilter.items()], key=lambda x:len(x[0]))
with open('alnDerep.tsv', 'w') as f:
    f.write('Refs\tCount\n')
    for item in alnSorted:
        f.write('{0}\t{1}\n'.format('|'.join(item[0]), item[1]))
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
content = []
for key, value in dist.items():
    content.append((key, value))
content = sorted(content, key=lambda x:x[0])
print(content)
with open('queryCount.tsv', 'w') as f:
    f.write('NumberAln\tCount\n')
    for line in content:
        f.write('{0}\t{1}\n'.format(line[0], line[1]))

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
    currentSubgraph = G.subgraph(nodes_to_remove).copy()
    sg.append((currentSubgraph, currentSubgraph.number_of_nodes()))
    G.remove_nodes_from(nodes_to_remove)
sg = sorted(sg, key=lambda x:x[1], reverse=True)
print('Graph divided into {0} subgraphs.'.format(len(sg)))
#%%
for item in sg:
    refNodes = nx.get_node_attributes(item, 'ref').keys()
    print(refNodes)