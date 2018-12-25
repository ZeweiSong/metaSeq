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
# Read in multiple alignments (usually two, but more than two is allowed) from the alignment files,
# and save normalized data based on the minimum depth
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
    # Node types: 1.reference, 2.forward query, 3.reverse query
    # Node atributes:
    #   attribute: 'r' (reference), 0 (forward), 1 (reverse)
    #   abundance: a tuple for query abundance (length = number of targets, usually = 2)
# Build the initial graph from the edge list (which is the alignment result)
# Add ref node list as attribute of the graph as a dictionary
def initGraph(alnNormalized):
    G = nx.Graph()
    refDict = {}
    for index, target in enumerate(alnNormalized):
        for alignment in target:
            ref = alignment[1]
            query = alignment[0]
            refDict[ref] = {}
            G.add_node(ref, attribute = 'r') # Assign attribute to a reference node
            G.add_node(query, attribute = index) # Assign attribute to a query node
            G.add_edge(ref, query) # Add the edge
    G.graph['ref'] = refDict # Graph attribute for the reference names
    G.graph['profile'] = {} # Graph attribute for the minimum profile
    G.graph['targetNumber'] = len(alnNormalized) # Graph attribute for the target number
    G = addAbundance(G, len(alnNormalized)) # Calculate the abundance for all reference node
    return G

# Add abundance value to the graph
    # Abundance value need to be updated after removal of the winner
def addAbundance(graph, n):
    for ref in graph.graph['ref'].keys():
        queryCount = [0] * n
        for query in graph.neighbors(ref): # Iterate all ref's neighbors
            queryCount[graph.nodes[query]['attribute']] += 1
        graph.nodes[ref]['abundance'] = tuple(queryCount)
    return graph

# Calculate the Effective Count value for current count string [m, n, ...]
    # This value can be considered as a weight for each reference.
def effectiveCount(countString):
    if len(countString) > 1:
        ave = sum(countString) / len(countString)
        stdev = statistics.stdev(countString)
        ec = ave - stdev
        if ec < 0:
            ec = 0.0
        return ec
    else:
        return countString[0]

# Return the new Graph after market competition
def competition(graph, greedy=True):
    import random
    # Get all references with ES >= 1
    references = graph.graph['ref'] # Dictionary for all available references
    refSurvivors = {} # Dictionary for saving survivors
    losers = [] # List for losers (ref with abundance < 1)
    for ref in references.keys():
        effScore = effectiveCount(graph.nodes[ref]['abundance'])
        if effScore >= 1.0:
            refSurvivors[ref] = effScore
        else:
            losers.append(ref)

    # Remove all losers from the graph
    for ref in losers:
        queries = graph.neighbors(ref)
        if not greedy:
            pool = [[]] * graph.graph['targetNumber']
            for q in list(queries):
                pool[graph.nodes[q]['attribute']].append(q) # allocate queries to corresponding target by their number
            minLength = min([len(i) for i in pool]) # Get the min size of the target
            querySurvived = []
            for target in pool:
                if len(target) > minLength:
                    random.shuffle(target)
                    target.sort(key=lambda x:graph.degree(x)[0])
                    querySurvived += target[minLength:]
            queries = [x for x in list(queries) if x not in querySurvived]
        else:
            pass
        graph.remove_nodes_from(list(queries))
        graph.remove_node(ref)
        del graph.graph['ref'][ref]

    # Check if the graph is empty
    if graph.number_of_edges() == 0:
        return graph
    else:
        # Get the winner (ref with largest EF)
        winner = sorted([i for i in refSurvivors.items()], key=lambda x:x[1], reverse=True)[0]

        # Remove winner from the graph
        ref = winner[0]
        queries = graph.neighbors(ref)
        graph.remove_nodes_from(list(queries)) # Remove all winner's query nodes
        graph.remove_node(ref) # Remove the winner node
        del graph.graph['ref'][ref]

        # Update the profile with the latest winner
        graph.graph['profile'][winner[0]] = winner[1]

        # Update the abundance string for the rest of references
        graph = addAbundance(graph, graph.graph['targetNumber'])

        return graph

#%%
# Define the winner take all function
def winnerTakeAll(aln, progress=False):
    G = initGraph(aln)
    if progress:
        print('Initial graph contains {0} references,'.format(len(G.graph['ref'])))
    while G.number_of_edges() > 0: # Need to change the condition to number of edge == 0
        G = competition(G)
        if progress:
            print('\t{0} references left, {1} are in the profile,'.format(len(G.graph['ref']), len(G.graph['profile'])))
    return G.graph['profile']

#%% Old logic for winner take all, save just for case
'''
# Pick the winner from the graph (with abundance tuple)
# TODO after the sorting, remove the zero abundance tail to speed up
def pickWinner(graph):
    # Iterate over all reference nodes:
    abundanceDict = {}
    for ref in graph.graph['ref'].keys():
        abundanceDict[ref] = effectiveCount(graph.nodes[ref]['abundance'])
    candidates = [i for i in abundanceDict.items() if i[1] >= 1] # Keep references with EF >= 1
    candidates.sort(key=lambda x:x[1], reverse = True)
    #print(candidates)
    return candidates
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
    losers = [i for i in graph.graph['ref'].items() if i[0] not in [j[0] for j in candidates]] #TODO need to fix this one
    print(losers)
    for item in losers:
        removeQueries = list(graph.neighbors(item[0]))
        graph.remove_nodes_from(removeQueries)
        graph.remove_node(item[0])
        del graph.graph['ref'][item[0]]
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
            print(candidates[0])
            G = removeWinner(G, candidates)
            G = addAbundance(G, targetNumber)
    profile = {}
    for ref in G.graph['profile'].keys():
        if G.graph['profile'][ref]: # This ref is not empty (broke)
            profile[ref] = G.graph['profile'][ref]
    return profile
'''