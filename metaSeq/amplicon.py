#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 09:51:54 2018

@author: Zewei Song
@email: songzewei@genomics.cn
"""
from __future__ import print_function
from __future__ import division
import networkx as nx
from metaSeq import io as seqIO
import random
import statistics

# Put all alignment into lists
# Read in multiple alignments (usually two) from the alignment files,
# and save normalized based on the minimum depth
def initAlignment(alnString):
    aln = []
    for alnFile in alnString:
        alnParser = seqIO.alignment(alnFile)
        alnCurrent = []
        for item in alnParser:
            alnCurrent.append(item)
        aln.append(alnCurrent)

    # Find the minimum depth (number of alignments)
    alnDepth = [len(i) for i in aln]
    minDepth = min(alnDepth)

    # Randomly sample all alignments into the minimum depth
    alnNormalized = []
    for item in aln:
        if len(item) > minDepth:
            alnNormalized.append(random.sample(item, minDepth))
        else:
            alnNormalized.append(item)
    return alnNormalized

#%% Build the alignment graph
    # Node: reference, forward reads, reverse reads
    # Node
    #   attribute: 'r' (reference), 0 (forward), 1 (reverse)
    #   abundance: a tuple for query abundance

# Build the initial graph from the edge list
# Add ref node list as attribute of the graph as a list
def initGraph(alnNormalized):
    G = nx.Graph()
    refDict = {}
    for index, target in enumerate(alnNormalized):
        for alignment in target:
            ref = alignment[1]
            query = alignment[0]
            refDict[ref] = {}
            G.add_node(ref, attribute = 'r')
            G.add_node(query, attribute = index)
            G.add_edge(ref, query)
    G.graph['ref'] = refDict
    G.graph['profile'] = refDict.copy()
    return G

# Add abundance value to the graph
def addAbundance(graph, n):
    for ref in graph.graph['ref'].keys():
        queryCount = [0] * n
        for query in graph.neighbors(ref): # Iterate all ref's neighbors
            queryCount[graph.nodes[query]['attribute']] += 1
        graph.nodes[ref]['abundance'] = tuple(queryCount)
    return graph

#%%
# Pick the winner from the graph (with abundance tuple)
def pickWinner(graph):
    # Iterate over all reference nodes:
    abundanceDict = {}
    for ref in graph.graph['ref'].keys():
        abundanceDict[ref] = effectiveCount(graph.nodes[ref]['abundance'])
    candidates = [i for i in abundanceDict.items()]
    candidates.sort(key=lambda x:x[1], reverse = True)
    return candidates[0]

#
# Calculate the Effective Count value for current count string [m, n, ...]
def effectiveCount(countString):
    ave = sum(countString) / len(countString)
    stdev = statistics.stdev(countString)
    ec = ave - stdev
    if ec < 0:
        ec = 0.0
    return ec

#%%
# Remove winner from the graph, remove all its neighbors
# Winner is the name of the reference
def removeWinner(graph, winner):
    removeQueries = list(graph.neighbors(winner))
    graph.remove_nodes_from(removeQueries)
    graph.remove_node(winner)
    del graph.graph['ref'][winner]
    return graph

# Return the winner take all profile
def winnerTakeAll(aln):
    G = initGraph(aln)
    targetNumber = len(aln)
    G = addAbundance(G, targetNumber)
    while not nx.is_empty(G):
        winner = pickWinner(G)
        G.graph['profile'][winner[0]] = winner[1]
        #print(G.graph['profile'])
        G = removeWinner(G, winner[0])
        G = addAbundance(G, targetNumber)
    profile = {}
    for ref in G.graph['profile'].keys():
        if G.graph['profile'][ref]: # This the ref is not empty (broke)
            profile[ref] = G.graph['profile'][ref]
    return profile