#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 16:59:32 2017

This script contains functions that take a sequence file as input.

@author: Zewei Song
@email: songzewei@genomics.cn
"""
from __future__ import print_function
from __future__ import division

# Function for checking the file type (gz, fasta, fastq)
# Return a tuple with two boolean value (T/F, T/F) --> (gz/not_gz, FASTQ/FASTA)
# FBI WARNING, this is a weak checker
    # It only check the file extension for gz file
    # It only check the first character for sequence type
def showMeTheType(filePath):
    import gzip
    seqType = {'>':False, '@':True}
    fileType = []
    if filePath.endswith('.gz'):
        fileType.append(True)
    else:
        fileType.append(False)
    if fileType[0]:
        with gzip.open(filePath, 'rt') as f:
            headSymbol = f.readline()[0]
    else:
        with open(filePath, 'r') as f:
            headSymbol = f.readline()[0]
    fileType.append(seqType[headSymbol])
    return tuple(fileType)


# an iterator oobject for reading a single sequence file. This is the most common Class.
# It will NOT check the format of the file, either can it deal with multiple line FASTA file.
# At my desktop computer, 1 M reads can be read in in about 2 seconds.
class sequence(object):
    def __init__(self, filePath):
        fileType = showMeTheType(filePath)
        if fileType[0]:
            import gzip
            self.file = gzip.open(filePath, 'rt')
        else:
            self.file = open(filePath, 'r')

        if fileType[1]:
            self.n = 4
        else:
            self.n = 2

    def __iter__(self):
        return self

    def __next__(self):
        record = []
        for i in range(self.n):
            line = self.file.readline().strip('\n')
            if line:
                record.append(line)
            else:
                raise StopIteration
        record[0] = record[0][1:]
        return record


# Same iterator but read in multiple record into memory at once
# By testing, read in file in trunk does not boost the reading speed.
class sequence_trunk(object):
    def __init__(self, filePath, trunk_size=100000):
        fileType = showMeTheType(filePath)
        if fileType[0]:
            import gzip
            self.file = gzip.open(filePath, 'rt')
        else:
            self.file = open(filePath, 'r')

        if fileType[1]:
            self.n = 4
        else:
            self.n = 2

    def __iter__(self):
        return self

    def __next__(self):
        record_trunk = []
        for record in range(self.trunk_size):
            record = []
            for i in range(self.n):
                line = self.file.readline().strip('\n')
                if line:
                    record.append(line)
                else:
                    if len(record_trunk) > 0:
                        return record_trunk
                    else:
                        raise StopIteration
            record[0] = record[0][1:]
            record_trunk.append(record)
        return record_trunk


# Iterator for two files
# It only work for files with ABSOLUTELY corresponding record.
class sequence_twin(object):
    def __init__(self, file_r1, file_r2):
        fileType = showMeTheType(file_r1)
        fileType_2 = showMeTheType(file_r2)
        if fileType[0] != fileType_2[0] or fileType[1] != fileType_2[1]:
            print('Inconsistent file type, are you serious?')
            self.n = 1
        if fileType[0]:
            import gzip
            self.r1 = gzip.open(file_r1, 'rt')
            self.r2 = gzip.open(file_r2, 'rt')
        else:
            self.r1 = open(file_r1, 'r')
            self.r2 = open(file_r2, 'r')
        if fileType[1]: self.n = 4
        else: self.n = 2

    def __iter__(self):
        return self

    def __next__(self):
        record = [[],[]]
        for i in range(self.n):
            line_r1 = self.r1.readline().strip('\n')
            line_r2 = self.r2.readline().strip('\n')
            if line_r1:
                record[0].append(line_r1)
                record[1].append(line_r2)
            else:
                raise StopIteration
        record[0][0] = record[0][0][1:]
        record[1][0] = record[1][0][1:]
        return record[0], record[1]


class sequence_twin_trunk(object):
    def __init__(self, file_r1, file_r2, trunk_size=100000):
        fileType = showMeTheType(file_r1)
        fileType_2 = showMeTheType(file_r2)
        if fileType[0] != fileType_2[0] or fileType[1] != fileType_2[1]:
            print('Inconsistent file type, are you serious?')
            self.n = 1
        if fileType[0]:
            import gzip
            self.r1 = gzip.open(file_r1, 'rt')
            self.r2 = gzip.open(file_r2, 'rt')
        else:
            self.r1 = open(file_r1, 'r')
            self.r2 = open(file_r2, 'r')
        if fileType[1]: self.n = 4
        else: self.n = 2
        self.trunk_size = trunk_size

    def __iter__(self):
        return self

    def __next__(self):
        r1_trunk = []
        r2_trunk = []
        for record in range(self.trunk_size):
            r1 = []
            r2 = []
            for i in range(self.n):
                line_r1 = self.r1.readline().strip('\n')
                line_r2 = self.r2.readline().strip('\n')
                if line_r1:
                    r1.append(line_r1)
                    r2.append(line_r2)
                else:
                    if len(r1_trunk) > 0:
                        return r1_trunk, r2_trunk
                    else:
                        raise StopIteration
            r1[0] = r1[0][1:]
            r2[0] = r2[0][1:]
            r1_trunk.append(r1)
            r2_trunk.append(r2)
        return r1_trunk, r2_trunk


# This function is still under test, it reads in file bytes, should be a bit faster than read in by line
# Need to find a way to remove header symbol
# This is 40% faster than reading by line.
class sequence_bytes(object):
    def __init__(self, filePath, size = 1000000,fastx='a'):
        self.file = open(filePath, 'r')
        self.size = size
        self.fastx = fastx
        self.tail = ''

        if fastx == 'a':
            self.n = 2
        elif fastx == 'q':
            self.n = 4
        else:
            print('Please specify the right format, "a" for FASTA and "q" for FASTQ.')
    def __iter__(self):
        return self

    def __next__(self):
        content = self.file.read(self.size)

        if content:
            content = self.tail + content
            self.tail = ''
            content = content.split('\n')
            self.tail = content[-1]
            tail_n = (len(content) - 1) % 4 # Set the line step to 4, set to 2 for FASTA
            if tail_n > 0:
                self.tail = '\n'.join(content[:-1][-1*tail_n::]) + '\n' + self.tail
            else:
                self.tail = self.tail
            #print(self.tail)
            content = content[:(len(content) - tail_n - 1)]
            content = [content[x:x+self.n] for x in range(0, len(content), self.n)]
            return content
        else:
            if self.tail:
                content = self.tail.split('\n')
                content = [content[x:x+self.n] for x in range(0, len(content), self.n)]
                return content
            else:
                raise StopIteration


# Write the content to a fastx file
def write_seqs(seq_content, filePath, fastx='a', mode='w'):
    count = 0
    if fastx == 'a':
        header = '>'
    elif fastx == 'q':
        header = '@'
    else:
        header = '-_-b'

    f = open(filePath, mode, newline='')
    for record in seq_content:
        label = header + record[0]
        for line in [label] + list(record[1:]):
            f.write('%s\n' % line)
        count += 1
    f.close()
    return count


#%% Functions for alignment
class alignment(object):
    def __init__(self, alnFile):
        self.fileName = alnFile
        self.aln = open(alnFile, 'r')

    def __iter__(self):
        return self

    def __next__(self):
        line = self.aln.readline()
        if line:
            line = line.strip('\n').split('\t')
            return line
        else:
            raise StopIteration

def readAlignment(alnFile, sort=True):
    aln = alignment(alnFile)
    alnList = []
    for line in aln:
        alnList.append(line)
    if sort: # Sort using the first column, query by default.
        alnList.sort(key=lambda x:x[0])
    return alnList

def writeAlignment(alnList, outFile):
    i = 0
    with open(outFile, 'w') as f:
        for line in alnList:
            f.write('{0}\n'.format('\t'.join([str(i) for i in line])))
            i += 1
    return i


#%%# Parse a sorted stLFR FASTA data by bead (barcode)
# Return a tuple containing all the short sequences in the bead
# Barcode is saved after the last "-", for example /102_1324_573
# Currently only support FASTA since all QC should be at the upstream
# TODO: rename this class to bead_fastx
# TODO: Parsing this way is too slow, need to think of alternatives
class stlfr_bead(object):
    def __init__(self, filePath, fastx='a'):
        if fastx == 'a':
            self.n = 2
        else:
            if fastx == 'q':
                self.n = 4
            else:
                print('Please specify the right FASTX format, a or q.')
                return None
        self.stop = False
        with open(filePath, 'r') as f:
            line1 = f.readline() # Read in the very first barcode
            self.barcode = line1.split('-')[0]
        self.file = open(filePath, 'r')
        self.current_bead = []

    def __iter__(self):
        return self

    def __next__(self):
        close = False
        while close == False:
            current_bead = self.current_bead
            record = []
            for i in range(self.n): # Read in the next sequence
                line = self.file.readline().strip('\n')
                if line:
                    record.append(line)
                else:
                    if self.stop:
                        raise StopIteration
                    else:
                        self.stop = True
                        current_bead = (self.barcode, tuple((tuple(i) for i in current_bead)))
                        return current_bead
            barcode = record[0].split('-')[0]
            if barcode == self.barcode:
                current_bead.append(record)
            else:
                current_bead = (self.barcode, tuple((tuple(i) for i in current_bead)))
                self.barcode = barcode
                self.current_bead = [record]
                close = True
                return current_bead


# Convert a stLFR sequence data (FASTQ) into JSON format per bead,
# Work for splited file with R1 and R2 files
# A dictionary is returned with barcode as key:
def fastq2json(fwdFile, revFile, barcodePosition=-2):
    beadJson = {}
    for r1, r2 in sequence_twin(fwdFile, revFile):
        barcode = r1[0].split('/')[barcodePosition]
        try:
            beadJson[barcode].append(r1 + r2)
        except KeyError:
            beadJson[barcode] = [r1+r2]
    return beadJson


# Convert the pair end merged files (one assembled, two unassembled) into bead dictionary
def mergepairs2bead(assemFile, fwdFile, revFile):
    bead = {} # Dictionary for storing bead
    ft = showMeTheType(assemFile)
    for record in sequence(assemFile):
        if ft[1]:
            barcode = record[0].split('/')[-2]
            try:
                bead[barcode].append([record[1], record[3]])
            except KeyError:
                bead[barcode] = [[record[1], record[3]]]
        else:
            barcode = record[0].split('/')[-2]
            try:
                bead[barcode].append([record[1]])
            except KeyError:
                bead[barcode] = [[record[1]]]
    ft1 = showMeTheType(fwdFile)
    ft2 = showMeTheType(revFile)
    if ft1 != ft2:
        print('Inconsisten forward and reverse file, come on, dude!')
    else:
        for r1, r2 in sequence_twin(fwdFile, revFile):
            if ft1[1]:
                barcode = r1[0].split('/')[-2]
                try:
                    bead[barcode].append([r1[1], r1[3], r2[1], r2[3]])
                except KeyError:
                    bead[barcode] = [[r1[1], r1[3], r2[1], r2[3]]]
            else:
                try:
                    bead[barcode].append([r1[1], r2[1]])
                except KeyError:
                    bead[barcode] = [[r1[1], r2[1]]]
        return bead


# Convert a pair of FASTQ record into list format as [seq1, q1, seq2, q2]
# All labels are discarded
def fastq2list(r1, r2):
    return [r1[1], r1[3], r2[1], r2[3]]


# Iterator for Bead from the JSON output
class bead_json(object):
    def __init__(self, filePath):
        import json
        self.file = open(filePath, 'r')
        self.head = []
        with open(filePath, 'r') as f:
            i = 1
            for line in f:
                if i <= 10:
                    i += 1
                    self.head.append(json.loads(line.strip('\n')))
    def __iter__(self):
        return self

    def __next__(self):
        import json
        line = self.file.readline().strip('\n')
        if line:
            return json.loads(line)
        else:
            raise StopIteration
#%%
# Return reverse compliment of a sequence
# This part is got from Stakoverflow
#(https://stackoverflow.com/questions/19570800/reverse-complement-dna) by corinna
# Works for Python 3
def revcomp(seq):
    return seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhdN', 'TGCAtgcaYRKMyrkmBVDHbvdhN'))[::-1]


# Count the line number of a file
# Count the line breaker \n by reading a segment in bytes
# The code is borrowed from Stackoverflow
def blocks(inputFile, size=65536):
    while True:
        b = inputFile.read(size)
        if not b: break
        yield b

def count_line(inputFile):
    with open(inputFile, "r", encoding="utf-8", errors='ignore') as f:
        return sum(bl.count("\n") for bl in blocks(f))