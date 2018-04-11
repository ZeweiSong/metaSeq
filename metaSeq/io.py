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
    def __init__(self, filePath, fastx='a', gzip='False', trunk_size=1):
        self.fastx = fastx
        self.gzip = gzip
        self.size = trunk_size
        if self.gzip:
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
        