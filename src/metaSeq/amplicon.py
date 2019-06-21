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
from itertools import combinations

# Reture a tuple of PCR barcode using for amplicon sequencing
# This list may be updated in the future for larger pool
def pcrBarcode():
    pcrbarcode = ('AACTGC', 'ACGGTT', 'AGCATG', 'AGGACT', 'CAGGAT', 'CAGTTC', 'CCATTG', \
                  'CGATAC', 'CTAGAG', 'CTCAGT', 'GAATCC', 'GACTTG', 'GATGAC', 'GCGTTA', \
                  'GTACGA', 'GTCCAT', 'TCCTAG', 'TCGATG', 'TGGCAT', 'TTGAGC', 'AACGAG')
    #pcrbarcode.sort()
    return pcrbarcode

# Reture the PCR barcide identified from the read(s)
# Reads is a tuple of read(s)
def getPCRBarcode(reads):
    if len(reads) == 1:
        barcode = reads[:6]
        return barcode
    else:
        bar1 = reads[0][:6]
        bar2 = reads[1][:6]
        if bar1 == bar2:
            return bar1
        else:
            return False


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
    alnDepth.sort()
    minDepth = min(alnDepth)

    print('Loaded {0} alignments/targets, containing {1} - {2} hits, ave = {3}'\
          .format(len(alnDepth), minDepth, alnDepth[-1], sum(alnDepth)/len(alnDepth)))

    # Randomly sample all alignments into the minimum depth
    alnNormalized = []
    for item in aln:
        if len(item) > minDepth:
            alnNormalized.append(random.sample(item, minDepth))
        else:
            alnNormalized.append(item)
    print('All targets normalized to the minimum depth: {0}.'.format(minDepth))
    return alnNormalized

# Count the number of reads in an alignment
def countAlignment(alnFile):
    aln = []
    alnParser = seqIO.alignment(alnFile)
    for item in alnParser:
        aln.append(item[0])
    aln = set(aln)
    return len(aln)
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
    #G['SRR3163732.100000']
    return G

# Build the reference graph, in which edge represents the number of overlapped queries
    # Now only work for single target alignment graph
def refGraph(alnGraph):
    G = nx.Graph()
    for node in alnGraph.nodes():
        if alnGraph.nodes[node]['attribute'] != 'r': # The node is not a ref, thus a query
            refs = list(alnGraph.neighbors(node))
            if len(refs) > 1:
                for pair in combinations(refs,2):
                    try:
                        G[pair[0]][pair[1]]['count'] += 1
                    except KeyError:
                        G.add_edge(pair[0],pair[1],count=1)
        else:
            G.add_node(node)
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
        return (ec, ave, stdev)
    else:
        return (countString[0], countString[0], 0)

# Alternatives to effectiveCount
# Median and Average maybe more suitable for technical replicates
def mediam(countString):
    if len(countString) > 1:
        ave = sum(countString) / len(countString)
        stdev = statistics.stdev(countString)
        return (statistics.median(countString), ave, stdev)
    else:
        return (countString[0], countString[0], 0)

def average(countString):
    if len(countString) > 1:
        ave = sum(countString) / len(countString)
        stdev = statistics.stdev(countString)
        return (ave, ave, stdev)
    else:
        return (countString[0], countString[0], 0)

# Return the new Graph after market competition
        # Consider Greedy and Less greedy methods
        # In Greedy mode, winner will grab all queries
        # In Less greedy mode, winner will balances the queries among targets to
        # maximize its EC, the extra queries are returned to other references
def competition(graph, greedy=True, weight='ec'):
    import random
    # Get all references with ES >= 1
    references = graph.graph['ref'] # Dictionary for all available references
    refSurvivors = {} # Dictionary for saving survivors
    losers = [] # List for losers (ref with abundance < 1)
    for ref in references.keys():
        if weight == 'ec': # Use Effective count
            effScore = effectiveCount(graph.nodes[ref]['abundance']) # need to add two other methods (medium and average)
        elif weight == 'median': # Use median number
            effScore = mediam(graph.nodes[ref]['abundance'])
        elif weight == 'average': # Use average number
            effScore = average(graph.nodes[ref]['abundance'])
        else:
            print('This is an new species of capitalist!')
        if effScore[0] >= 1.0:
            refSurvivors[ref] = effScore
        else:
            losers.append(ref)

    # Get the winner (ref with largest EF)
    winner = sorted([i for i in refSurvivors.items()], key=lambda x:x[1][0], reverse=True)
    if len(winner) > 0:
        winner = winner[0]
    else:
        graph.remove_nodes_from(references) # There are nodes left, but none is eligible winner candidate
        return graph

    winner = list(winner)
    winner[1] = list(winner[1])
    print('{0} is the winner (weight on {3} = {4}), with ave = {1}, stdev = {2}.'.format(winner[0], winner[1][1], winner[1][2], weight, winner[1][0]))

    # Remove winner queires
        # Has two methods: Greedy and Less greedy
        # In the Less greedy mehtod, queries in winner are balanced, and extra queries are recycled back to the graph
    if greedy: # The greedy method, by default
        # Remove winner from the graph
        ref = winner[0]
        queries = graph.neighbors(ref)
        graph.remove_nodes_from(list(queries)) # Remove all winner's query nodes
        graph.remove_node(ref) # Remove the winner node
        del graph.graph['ref'][ref]
        # CRITICAL: need to update ref abundance string after each round
        # Update the abundance string for the rest of references
        graph = addAbundance(graph, graph.graph['targetNumber'])
    else: # The less greedy method, queries are recycled to maximize the balance of winner
            # for nTarget = 2, keep the lower depth
            # for nTarget > 2, need to think a different way (use medium depth, average depth ...)
        ref = winner[0]
        # allocate queries to corresponding target by their number
        queries = graph.neighbors(ref)
        pool = {}
        for i in range(graph.graph['targetNumber']):
            pool[i] = []
        for q in list(queries):
            pool[graph.nodes[q]['attribute']].append(q)
        minLength = min([len(i) for i in pool.values()]) # Get the min size of the target
        querySurvived = []
        for target in pool.values():
            if len(target) > minLength:
                random.shuffle(target) # shuffle the query list
                target.sort(key=lambda x:graph.degree(x)) # Sort query by their degree
                querySurvived += target[minLength:]
        print('\t{0} queries returned to the graph by the less greedy {1},'.format(len(querySurvived), winner[0]))
        queries = [x for x in list(queries) if x not in querySurvived]
        graph.remove_nodes_from(list(queries)) # Remove all winner's query nodes
        graph.remove_node(ref) # Remove the winner node
        del graph.graph['ref'][ref]
        winner[1][0] = minLength
        # Update the abundance string for the rest of references
        graph = addAbundance(graph, graph.graph['targetNumber'])

    # Remove all losers from the graph
    for ref in losers:
        if greedy:
            queries = graph.neighbors(ref)
            graph.remove_nodes_from(list(queries))
            graph.remove_node(ref)
            del graph.graph['ref'][ref]
        else: # Recheck the EC value of these losers, since some queries got recycled in the Less greedy method
            if effectiveCount(graph.nodes[ref]['abundance'])[0] < 1:
                queries = graph.neighbors(ref)
                graph.remove_nodes_from(list(queries))
                graph.remove_node(ref)
                del graph.graph['ref'][ref]
            else:
                pass

    # Update the profile with the latest winner
    graph.graph['profile'][winner[0]] = winner[1][0]
    return graph
#%%
# Define the winner take all function
def winnerTakeAll(aln, progress=True, greedy=False, weight='ec'):
    G = initGraph(aln)
    if progress:
        print('Initial graph contains {0} references,'.format(len(G.graph['ref'])))
    while G.number_of_edges() > 0: # Need to change the condition to number of edge == 0
        G = competition(G, greedy=greedy, weight=weight)
        if progress:
            print('\t{0} references left, {1} are in the profile,'.format(len(G.graph['ref']), len(G.graph['profile'])))
    return G.graph['profile']