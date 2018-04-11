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
    
    def __iter__(self):
        return self
    
    def __next__(self):
        record_trunk = []
        for record in range(self.trunk_size):
            record = []
            for i in range(self.n):
                line = self.file.readline().strip('\n')
                if line:
                    print(line)
                    record.append(line)
                    print(record)
                else:
                    if len(record_trunk) > 0:
                        return record_trunk
                    else:
                        raise StopIteration
            record[0] = record[0][1:]
            record_trunk.append(record)
        return record_trunk