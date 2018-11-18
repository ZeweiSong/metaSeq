#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 15:24:51 2018

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%%
from __future__ import print_function
from __future__ import division
import networkx as nx
from itertools import combinations
import random
import matplotlib.pyplot as plt
import community # Louvain Community Detection (https://github.com/taynaud/python-louvain/)
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', help='Input distance file')
args = parser.parse_args()

distFile = args.i

# Read distance data as edges
G=nx.Graph()
with open(distFile, 'r') as f:
    for line in f:
        line = line.strip('\n').split('\t')
        G.add_weighted_edges_from([(line[0],line[1],line[2])])

print('Graph contains {0} nodes and {1} edges.'.format(G.number_of_nodes(), G.number_of_edges()))
#nx.draw(G) # Can not draw in linux OS

# Find the modules in the graph
#%%Finds communities in a graph using the Girvan–Newman method
# This method generate modules (communities) in different levels
import networkx.algorithms.community as nxc
c = nxc.girvan_newman(G)
comPool = []
for item in c:
    comPool.append(tuple(sorted(item)))
print(len(comPool))
print(comPool[0])
print(comPool[49])
print(comPool[90])

#%% Returns communities in G as detected by Fluid Communities algorithm
c = nxc.asyn_fluidc(G, 4) # 4 is the number of expected communities
for item in c:
    print(item)

#%% Greedy modularity
c = list(nxc.greedy_modularity_communities(G))
print(c)

#%%
# Find k-clique communities in graph using the percolation method.
# Definition of cliques: https://en.wikipedia.org/wiki/Clique_(graph_theory)
c = list(nxc.k_clique_communities(G, 3))
print(c)

#%% Louvain Community Detection
partition = community.best_partition(G)

#%% Calculate the modularity of this graph
print('The modularity of G is: {0}.'.format(nxc.modularity(G, comPool[49])))
