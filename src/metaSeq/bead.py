#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 16:59:32 2017

@author: songz
"""
from __future__ import print_function
from __future__ import division
from metaSeq import qc as seqQC
from metaSeq import io as seqIO
import json


# Create an iterator Class that iterates on stLFR sequence file.
# Iterate by grouping sequences with same barcode.
class beadSequenceIterator(object):
    def __init__(self, sequenceFile):
        pass


# Create a bead Class
# Convert a JSON loaded dictionary into a Bead Class
# TODO need to add option for FASTA and FASTQ file
class beadSequence(object):
    def __init__(self, jsonDict):
        self.json = jsonDict
        self.barcode = list(jsonDict.keys())[0] # barcode associated with this bead
        self.fragments = list(jsonDict.values())[0] # fragments associated with this bead
        self.assembled = [] # assembled sequences (1 in each list)
        self.unassembled = []  # unassembled pair sequences (2 in each list)
        for record in self.fragments:
            if len(record) == 1:
                self.assembled.append(record)
            elif len(record) == 2:
                self.unassembled.append(record)

    # Generate a FASTA list
    # ((barcode/0, sequence), (barcode/1, sequence))
    # List can be wrote to FASTA file using metaseq.io.write_seqs()
    # /9/ for assembled read
    # /1/ and /2/ for unassembled paired reads
    def fastaSequences(self):
        fasta = []
        i = 0
        for record in self.assembled:
            fasta.append((self.barcode + '/9/' + str(i), record[0]))
            i += 1
        for record in self.unassembled:
            fasta.append((self.barcode + '/1/' + str(i), record[0]))
            fasta.append((self.barcode + '/2/' + str(i), record[1]))
            i += 1
        return fasta


    # Write the sequecnes into fasta file with barcode as file name
    def fastaWrite(self, folder=''):
        fasta = self.fastaSequences()
        fileName = folder + self.barcode + '.fa'
        count = seqIO.write_seqs(fasta, fileName)
        return count


    # Write the bead into a JSON file
    def jsonWrite(self, outputFile, mode='a'):
        with open(outputFile, mode, newline='') as f:
            f.write('%s\n' % json.dumps({self.barcode: self.assembled + self.unassembled}))



# A class for kmer group of a given bead
class beadKmer(object):
    def __init__(self, beadSequence):
        self.sequence = []
        pass

    # Reture the kmer set of this bead
    def kmer(self, kmer_size=21):
        pass


class kmerDistance(object):
    def __init__(self, kmer1, kmer2):
        self.kmer1 = kmer1
        self.kmer2 = kmer2

    def jaccard(self):
        pass

    def distance(self):
        pass

#%% Create an iterator that iterate over alignment output of stLFR data
# Iterate by grouping alignments with same barcode
#  I recommend to use this format based on Vsearch: query+target+ql+tl+id+qilo+qihi+tilo+tihi
# WARNING output of Vsearch changed the order of alignment! You need to re-order the alignment before iterate
class beadAlignmentIterator(object):
    def __init__(self, alnFile):
        with open(alnFile, 'r') as f:
            line1 = f.readline() # Read in the very first barcode
            self.barcode = line1.strip('\n').split('\t')[0].split('/')[0]
        self.stop = False
        self.current_bead = []
        self.aln = open(alnFile, 'r')

    def __iter__(self):
        return self

    def __next__(self):
        close = False
        while close == False:
            current_bead = self.current_bead
            record = []
            line = self.aln.readline().strip('\n')
            if line:
                record.append(line.split('\t'))
            else:
                if self.stop:
                    raise StopIteration
                else:
                    self.stop = True
                    current_bead = (self.barcode, tuple((tuple(i) for i in current_bead)))
                    return current_bead
            barcode = line.strip('\n').split('\t')[0].split('/')[0]
            if barcode == self.barcode:
                current_bead.append(line.split('\t'))
            else:
                current_bead = (self.barcode, tuple((tuple(i) for i in current_bead)))
                self.barcode = barcode
                self.current_bead = [line.split('\t')]
                close = True
                return current_bead


# An Object for basic bead based alignment information
class beadAlignment(object):
    def __init__(self, beadAln):
        self.barcode = beadAln[0]
        self.alnRaw = beadAln[1]
        self.fragmentCount = {}

        # For the input alignemnt, if only one read from a pair alignment to a reference, we need to remove it.
        alnDict = {} # create a dict with ref and qyuery label as key
        for line in self.alnRaw:
            query = line[0]
            ref = line[1]
            try:
                alnDict[ref][query] = line
            except KeyError:
                alnDict[ref] = {query: line}

        # Find the pair to be removed (paired read is label as /1/# and /2/#)
        keepDict = {} # Store all sequence lables to keep
        for key, value in alnDict.items(): # Find unpaired alignment reads per reference
            fragmentDict = {} # Dict to coun the number of read in each pair
            keepDict[key] = [] # Add current reference to the keep Dict
            queryList = list(value.keys()) # All query names of current reference
            self.fragmentCount[key] = 0 # Count the fragment for current reference
            for item in queryList: # Parse each sequence
                fragmentSplit = item.split('/')
                fragmentType = fragmentSplit[1]
                fragmentNumber = fragmentSplit[2]
                if fragmentType == '9': # This is an assembled read
                    keepDict[key].append(item)
                    self.fragmentCount[key] += 1
                else: # This is a paired read
                    try:
                        fragmentDict[fragmentNumber] += 1
                    except KeyError:
                        fragmentDict[fragmentNumber] = 1
            for number, count in fragmentDict.items():
                if count == 2: # Both read in a pair have eligible alignment in current reference
                    keepDict[key].append(self.barcode + '/1/' + number)
                    keepDict[key].append(self.barcode + '/2/' + number)
                    self.fragmentCount[key] += 1
                else:
                    pass

        # Store all eligible read into a new list
        self.aln = []
        for key, value in keepDict.items():
            for item in value:
                self.aln.append(alnDict[key][item])

    def references(self):
        refDict = {}
        for line in self.aln:
            try:
                refDict[line[1]] += 1
            except KeyError:
                refDict[line[1]] = 1
        return refDict

    def queries(self):
        queryDict = {}
        for line in self.aln:
            try:
                queryDict[line[0]] += 1
            except KeyError:
                queryDict[line[0]] = 1
        return queryDict

    def minSet(self):
        return winnerTakeAll(self.aln)
#%%

# An Class object for the minimum set alignment (output of beadAlignment.minSet())
# Alignment is in the format of query+target+ql+tl+id+qilo+qihi+tilo+tihi (Vsearch userfields)
class beadMinSet(object):
    def __init__(self, minSet):
        self.minSet = minSet
        self.barcode = list(list(self.minSet.values())[0].keys())[0].split('/')[0]
        self.references = list(self.minSet.keys()) # List of all reference associated with this bead
        self.refLength = {} # dictionary for reference and its length in bp
        for key, value in self.minSet.items(): # Length of each reference, Count of fragment (need to fix) of each reference
            queries = list(value.keys())
            firstQuery = queries[0]
            refL = int(value[firstQuery][1])
            self.refLength[key] = refL

    # Count the number of fragment (NOT read) per reference
    def fragmentCount(self):
        refCount = {} # A dictionary for reference and number of alignments it contains
        for key, value in self.minSet.items(): # Length of each reference, Count of fragment (need to fix) of each reference
            queries = list(value.keys())
            fragmentList = []
            for item in queries:
                item = item.split('/')[1]
                if item in ['9','1']:
                    fragmentList.append(item)
            refCount[key] = len(fragmentList) # Need to consider paired read (count as one)
        return refCount

    # Return all position number using an alignment string
    def positionString(self, aln):
        start = int(aln[5])
        end = int(aln[6])
        posList = []
        for i in range(end - start + 1):
            posList.append(start + i)
        return posList

    # Return the covered position of a given reference
    def refCoveredBases(self, refName):
        alnList = list(self.minSet[refName].values())
        basePosition = []
        for aln in alnList:
            basePosition += self.positionString(aln)
        return list(set(basePosition))

    # A report on the coverage of eacch reference
    def refCoverageReport(self):
        report = {}
        for ref in self.references:
            coveredBase = len(self.refCoveredBases(ref))
            refL = self.refLength[ref]
            report[ref] = [coveredBase, refL, coveredBase / refL]
        return report

#%%
# Remove duplicated reads from a single bead Class
# The reverse complimentary sequence is considered
# Return a new dictionary (used as input of beadSequence())
def derep(bead):
    seqs = list(bead.values())[0]
    derep_dict = {}
    for record in seqs:
        if len(record) == 1: # assembled read
            key = [(record[0],), (seqIO.revcomp(record[0]),)]
            key.sort(key=lambda x:x[0])
            key = tuple(key)
            #print(key)
            try:
                derep_dict[key] += 1
            except KeyError:
                derep_dict[key] = 1
        if len(record) == 2: # unassembled reads
            key = [(record[0], record[1]), (seqIO.revcomp(record[1]), seqIO.revcomp(record[0]))]
            key.sort(key=lambda x:x[0])
            key = tuple(key)
            try:
                derep_dict[key] += 1
            except KeyError:
                derep_dict[key] = 1
    barcode = list(bead.keys())[0]
    bead_dict = {}
    for key, value in derep_dict.items():
        bead_dict[key[0]] = value
    return {barcode: bead_dict}


# Remove low quality reads from a single bead Class
# Return a new bead Class
# Now always return without quality score (has a return_format value for the future)
def maxEE(bead, maxee=1, return_format='a'):
    seqs = list(bead.values())[0]
    seqs_qc = []
    p_dict = seqQC.qual_score()
    for record in seqs:
        if len(record) == 2: # assembled read
            if seqQC.ee(record[1], p_dict) < maxee:
                seqs_qc.append([record[0]])
        elif len(record) == 4: # unassembled read1 and read2
            if seqQC.ee(record[1], p_dict) < maxee and seqQC.ee(record[3], p_dict) < maxee:
                seqs_qc.append([record[0], record[2]])
    barcode = list(bead.keys())[0]
    return {barcode: seqs_qc}


#%% Winner take all (wta) methods for finding the minimum set of alignments (all possible hits)
# Input is alignment results saved in tuple as (query, target, ... , tilo, tihi)
# Return a dict use reference as key, and all associated alignment as value
def winnerTakeAll(aln):
    minSet = {}
    alnDict = createDict(aln)
    i = 0
    #print(rdict)
    while len(alnDict) > 0:
        i += 1
        #print('Round {0}'.format(i))
        #print(len(rdict))
        winnerAln = pickWinner(alnDict)
        #print(winner)
        winner = winnerAln[0]
        minSet[winner] = alnDict[winner]
        alnDict = removeWinner(alnDict, winnerAln)
    return minSet


# Create a ref-query dictionary
def createDict(aln):
    alnDict = {}
    for line in aln:
        ref = line[1]
        query = line[0]
        try:
            alnDict[ref][query] = tuple(line[2:])
        except KeyError:
            alnDict[ref] = {query: tuple(line[2:])}
    return alnDict

# From the current alignment pick a winner (ref with the most query)
def pickWinner(alnDict):
    candidates = [[key, len(value)] for key, value in alnDict.items()]
    candidates.sort(key=lambda x:x[1], reverse=True)
    winner = candidates[0][0]
    aln = []
    for key, value in alnDict[winner].items():
        aln.append((key,) + value)
    aln = tuple(aln)
    return (winner, aln)

# Remove the winner from the alignment dictionary
def removeWinner(d, winnerAln):
    winner = winnerAln[0]
    #trash = d.pop(winner, None)
    trash = d[winner]
    del d[winner]
    for key in d.keys():
        for query in trash.keys():
            #a = d[key].pop(query, None)
            try:
                del d[key][query]
            except KeyError:
                pass
    losers = []
    for key in d.keys():
        if len(d[key]) == 0: # Current reference is empty
            losers.append(key)
    for loser in losers:
        #a = d.pop(loser, None)
        del d[loser]
    return d