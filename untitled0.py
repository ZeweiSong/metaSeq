#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 15:27:36 2018

Change the unite sequence file label to SH only, and save the taxa in a single file

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%%
from __future__ import print_function
from __future__ import division
from metaSeq import io as seqIO

handle = seqIO.sequence('ee_its_database.fasta', fastx='a')
unite = []
for item in handle:
    unite.append(item)

fasta = []
taxa = []
for item in unite:
    label = item[0].split('|')
    newLabel = label[0] + '|' + label[2]
    fasta.append((newLabel, item[1]))
    taxa.append((newLabel, label[4]))

seqIO.write_seqs(fasta, 'ee_its_sequences.fa', fastx='a')
with open('ee_its_taxonomy.txt', 'w') as f:
    for line in taxa:
        line = '\t'.join(line)
        f.write('{0}\n'.format(line))