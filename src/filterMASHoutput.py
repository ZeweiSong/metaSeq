#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 15:10:37 2018

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%%
from __future__ import print_function
from __future__ import division
inputFile = 'dist.txt'
pair = []
threshold = 0.04
with open(inputFile, 'r') as f:
    for line in f:
        line = line.strip('\n').split('\t')
        source = line[0].split('/')[1]
        target = line[1].split('/')[1]
        distance = float(line[2])
        
        if source != target and 0.02 <= distance <= threshold:
            pair.append([source, target, distance])
with open('dist_002_004.txt', 'w') as f:
    for line in pair:
        f.write('{0}\n'.format('\t'.join([str(i) for i in line])))