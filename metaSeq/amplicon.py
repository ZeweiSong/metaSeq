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
    G.graph['profile'] = {}
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
# TODO after the sorting, remove the zero abundance tail to speed up
def pickWinner(graph):
    # Iterate over all reference nodes:
    abundanceDict = {}
    for ref in graph.graph['ref'].keys():
        abundanceDict[ref] = effectiveCount(graph.nodes[ref]['abundance'])
    candidates = [i for i in abundanceDict.items() if i[1] >= 1] # Keep references with EF >= 1
    candidates.sort(key=lambda x:x[1], reverse = True)
    return candidates

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
# TODO Also need to remove all references with EF < 1
def removeWinner(graph, candidates):
    # STEP1. Remove the winner
    removeQueries = list(graph.neighbors(candidates[0][0]))
    graph.remove_nodes_from(removeQueries)
    graph.remove_node(candidates[0][0])
    del graph.graph['ref'][candidates[0][0]]

    # STEP2. Remove the losers (references not in cadidates list)
    losers = [i for i in graph.graph['ref'].items() if i[0] not in candidates]
    for item in losers:
        removeQueries = list(graph.neighbors(item))
        graph.remove_nodes_from(removeQueries)
        graph.remove_node(item)
        del graph.graph['ref'][item]
    return graph

# Return the winner take all profile
# TODO need to change according to pickWinner()
def winnerTakeAll(aln):
    G = initGraph(aln)
    targetNumber = len(aln)
    G = addAbundance(G, targetNumber)
    while not nx.is_empty(G):
        candidates = pickWinner(G)
        if candidates[0][1] <= 0: # Stop if the winner is ZERO (or maybe 1?)
            break
        else:
            G.graph['profile'][candidates[0][0]] = candidates[0][1]
            #print(winner)
            G = removeWinner(G, candidates)
            G = addAbundance(G, targetNumber)
    profile = {}
    for ref in G.graph['profile'].keys():
        if G.graph['profile'][ref]: # This ref is not empty (broke)
            profile[ref] = G.graph['profile'][ref]
    return profile