#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 10:50:36 2018

@author: Zewei Song
@email: songzewei@genomics.cn
"""
from __future__ import print_function
from __future__ import division
from metaSeq import io as seqIO
import random
#%%
# Return the kmer list of a given sequences
# It will return all kmer in both strand WITHOUT dereplication
# So the kmer will be REDUNDANT in the output
def kmer(seq, size):
    seq = seq.upper()
    seq_length = len(seq)
    kmer = []
    for i in range(seq_length - size + 1):
        kmer.append(seq[i:i+size])
    seq_rev = seqIO.revcomp(seq)
    for i in range(seq_length - size + 1):
        kmer.append(seq_rev[i:i+size])
    return kmer

# Return a kmer table with abundance on the given kmer length for a set of sequences
# The input sequences is all in a list or tuple (seq1, seq2, seq3)
def kmerCount(seqs, size):
    kmer_list = [item for seq in seqs for item in kmer(seq, size)]
#    for item in seqs:
#        kmer_list += kmer(item, size)
    kmer_table = {}
    for item in kmer_list:
        kmer_table[item] = kmer_table.get(item, 1) + 1 # one line code for adding new keys
    return kmerTable(kmer_table) # return a kmerTable class


# A class for the kmer count of a given sequence set
class kmerTable(object):
    def __init__(self, kmerTable):
        self.kmers = list(kmerTable.keys())
        self.kmerSize = len(self.kmers[0])
        self.set = set(self.kmers)
    
    # Reture the minimizer of current kmer set (first kmer sorted)
    def minimizer(self):
        self.kmers.sort()
        return self.kmers[0]
    
    # Random sample the kmer pool
    def randomSample(self, sampleSize):
        if sampleSize <= len(self.kmers):
            return random.sample(self.kmers, sampleSize)
        else:
            return None


# Return distance of two kmer pool
# Input is a pair of lists/tuples, each contains a kmer pool
    # This is usually the output of itertools.combinations
class kmerDistance(object):
    def __init__(self, kmerPools):
        self.poolA = kmerPools[0]
        self.poolB = kmerPools[1]
    
    # return the jarccard index
    def jaccard(object):
        total_size = len(kmer_table)
        shared_size = 0
        for line in kmer_table:
            if line[1] > 0 and line[2] > 0: # Current kmer is present in both sets
                shared_size += 1
        return shared_size / total_size