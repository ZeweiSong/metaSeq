#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 16:59:32 2017

@author: songz
"""
from __future__ import print_function
from __future__ import division

# an iterator oobject for reading a single sequence file
# It will NOT check the format of the file, either can it deal with multiple line FASTA file.
class sequence(object):
    def __init__(self, filePath, fastx='a', gz=False, trunk_size=1):
        self.fastx = fastx
        self.gzip = gz
        self.size = trunk_size
        if self.gzip:
            import gzip
            self.file = gzip.open(filePath, 'rt')
        else:
            self.file = open(filePath, 'r')
        if fastx == 'a':
            self.n = 2
        elif fastx == 'q':
            self.n = 4
        else:
            print('Please specify the right format, "a" for FASTA and "q" for FASTQ.')
    
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
class sequence_trunk(object):
    def __init__(self, filePath, fastx='a', gz=False, trunk_size=2):
        self.fastx = fastx
        self.gzip = gz
        self.trunk_size = trunk_size
        if self.gzip:
            import gzip
            self.file = gzip.open(filePath, 'rt')
        else:
            self.file = open(filePath, 'r')
        if fastx == 'a':
            self.n = 2
        elif fastx == 'q':
            self.n = 4
        else:
            print('Please specify the right format, "a" for FASTA and "q" for FASTQ.')
            self.n = 1
    
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
    def __init__(self, file_r1, file_r2, fastx='a', gz=False):
        self.fastx = fastx
        self.gzip = gz
        if self.gzip:
            import gzip
            self.r1 = gzip.open(file_r1, 'rt')
            self.r2 = gzip.open(file_r2, 'rt')
        else:
            self.r1 = open(file_r1, 'r')
            self.r2 = open(file_r2, 'r')
        if fastx == 'a': self.n = 2
        elif fastx == 'q': self.n = 4
        else:
            print('Please specify the right format, "a" for FASTA and "q" for FASTQ.')
            self.n = 1
    
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
        return record


# Write the content to a fastq file
def write_seqs(seq_content, filePath, fastx='a', gz=False):
    count = 0
    if fastx == 'a':
        n = 2
        header = '>'
    elif fastx == 'q':
        n = 4
        header = '@'
    else:
        n = 1
        header = ''
    
    if not gz:
        f = open(filePath, 'w')
    else:
        import gzip
        f = gzip.open(filePath, 'w')
    for record in seq_content:
        label = header + record[0]                
        for line in [label] + record[1:]:
            f.write('%s\n' % line)
            count += 1
    f.close()
    return count