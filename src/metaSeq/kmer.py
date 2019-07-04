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
import numpy as np
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
def kmerCount(bead, size):
    seqs = [i[1] for i in bead.fastaSequences()]
    kmer_list = [item for seq in seqs for item in kmer(seq, size)]
#    for item in seqs:
#        kmer_list += kmer(item, size)
    kmer_table = {}
    for item in kmer_list:
        kmer_table[item] = kmer_table.get(item, 1) + 1 # one line code for adding new keys
    return kmerTable(kmer_table, bead.barcode) # return a kmerTable class


# A class for the kmer count of a given sequence set
class kmerTable(object):
    def __init__(self, kmer_table, barcode):
        self.kmers = list(kmer_table.keys())
        self.kmerSize = len(self.kmers[0])
        self.set = set(self.kmers)
        self.barcode = barcode
    
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
# Input is a pair of lists/tuples, from class.self.set, or read in from JSON kmer file
    # This is usually the output of itertools.combinations
class kmerDistance(object):
    # Usage
    # kd = kmerDistance([kmer1, kmer2])
    # md = kd.mashDistance
    def __init__(self, kmerPools):
        self.poolA = kmerPools[0]
        self.poolB = kmerPools[1]
        self.kmerLength = len(list(self.poolA)[0])
    
    # Return overlap of two kmerSet
    def overlap(self, smallSet, largeSet):
        overlap = []
        for item in largeSet:
            if item in smallSet:
                overlap.append(item)
        return tuple(overlap)
    
    # return the jarccard index
    def jaccard(self):
        lengthA = len(self.poolA)
        lengthB = len(self.poolB)
        if lengthA <= lengthB:
            overlapKmer = self.overlap(self.poolA, self.poolB)
        else:
            overlapKmer = self.overlap(self.poolB, self.poolA)
        lengthShared = len(overlapKmer)
        jaccardIndex = lengthShared / (lengthA + lengthB - lengthShared)
        return jaccardIndex
    
    # Return Mash distance of two kmer Sets
    def mashDistance(self):
        jac  = self.jaccard()
        if jac == 0:
            return 1
        else:
            mashD = (-1/self.kmerLength) * np.log((2 * jac)/(1 + jac))
            return mashD
    
    # Return the kmer distance, accounted for kmer abundance, usng Euclidean distance
    # May need to standardize before calculating
    def kmer_euclidean(set1, set2, size):
        pass