#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 16:59:32 2017

@author: songz
"""
#%%
from __future__ import print_function
from __future__ import division

class fastq_parser(object):
    # Define a iterator for reading the FASTQ file by record
    # The iterator will NOT check the format of input FASTQ file
    def __init__(self, filepath):
        self.file = open(filepath, 'r') # File object

    def __iter__(self):
        return self
    
    def __next__(self):
        record = []
        for i in range(4): # Read in four lines for iteration
            line = self.file.readline().strip('\n')
            if line:
                record.append(line)
            else:
                raise StopIteration
        record[0] = record[0][1:]
        return record

class fastq_parser_gz(object):
    # Define a iterator for reading the FASTQ file in Gzip format by record
    # The iterator will NOT check the format of input FASTQ file
    def __init__(self, filepath):
        import gzip
        self.file = gzip.open(filepath, 'r') # Gzip file object
        
    def __iter__(self):
        return self
    
    def __next__(self):
        record = []
        for i in range(4): # Read in four lines for iteration
            line = self.file.readline().decode('utf-8').strip('\n')
            if line:
                record.append(line)
            else:
                raise StopIteration
        record[0] = record[0][1:]
        return record

class read_bytes(object):
    # Read a file by trunk with n bytes
    def __init__(self, filepath, n):
        self.file = open(filepath, 'r')
        self.trunk_size = n
    
    def __iter__(self):
        return self
    
    def __next__(self):
        trunk = self.file.read(self.trunk_size)
        if trunk == '':
            raise StopIteration
        else:
            return trunk

class fastq_parser_bt(object):
    # Define an iterator that read in a FASTQ file by trunk of bytes
    def __init__(self, filepath, n):
        self.file = open(filepath, 'r')
        self.trunk_size = n
        self.tail = ''
    def __iter__(self):
        return self
    
    def __next__(self):
        trunk = self.file.read(self.trunk_size)
        if trunk == '': 
            raise StopIteration
        else:
            trunk = self.tail + trunk # Add the previous tail
            trunk = trunk.split('\n')
            tail_n = len(trunk) % 4
            self.tail = ''.join(trunk[-1*tail_n::]) # Need to deal with the \n problem
            return trunk

class fastq_double_gz(object):
    def __init__(self, file1, file2):
        import gzip
        self.file1 = gzip.open(file1, 'r')
        self.file2 = gzip.open(file2, 'r')
    
    def __iter__(self):
        return self
    
    def __next__(self):
        r1 = []
        r2 = []
        for i in range(4):
            r1_line = self.file1.readline().decode('utf-8').strip('\n')
            r2_line = self.file2.readline().decode('utf-8').strip('\n')
            if r1_line:
                r1.append(r1_line)
                r2.append(r2_line)
            else:
                raise StopIteration
        r1[0] = r1[0][1:]
        r2[0] = r2[0][1:]
        return [r1, r2]


# Write the content to a fastq file
def write_fastq(seq_content, filepath, gz = False):
    count = 0
    if not gz:
        with open(filepath, 'w') as f:
            for record in seq_content:
                label = '@' + record[0]                
                for line in [label] + record[1:]:
                    f.write('%s\n' % line)
                count += 1
        return count
    else:
        import gzip
        with gzip.open(filepath, 'w') as f:
            for record in seq_content:
                label = '@' + record[0]
                for line in [label] + record[1:]:
                    output_line = line + '\n'
                    f.write(output_line.encode())
                count += 1
        return count