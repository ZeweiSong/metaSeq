#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 10:36:47 2018

Coders who love to comment their code are unlikely to have bad luck.

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%%
from __future__ import print_function
from __future__ import division
import numpy as np
n = 1036**3
r = 50 * 10**6
print('Pool is {0:8.2f} M, draw is {1:8.2f} M.'.format(n/10**6, r/10**6))
print('Draw {0} beads from the pool.'.format(r))

for i in range(1000):
    draw = np.random.choice(n, r, replace=True)
    
    print('Counting identical beads.')
    count = {}
    for item in draw:
        count[item] = count.get(item, 0) + 1
    print(len(count))
    print('Getting the distribution on singleton and non-singleton beads.')
    dist = {}
    for key, value in count.items():
        dist[value] = dist.get(value, 0) + 1

    print('Writing {0} report'.format(i))
    dist = sorted([(i, j) for i,j in dist.items()], key=lambda x:x[0])
    with open('draw/beadDist_{0}.tsv'.format(i), 'w') as f:
        f.write('Repeats\tCount\n')
        for line in dist:
            f.write('{0}\t{1}\n'.format(line[0], line[1]))

#%%
dist = {1:[], 2:[], 3:[], 4:[], 5:[], 6:[]}
for i in range(100):
    d = {1:0, 2:0, 3:0, 4:0, 5:0, 6:0}
    with open('draw/beadDist_{0}.tsv'.format(i), 'r') as f:
        f.readline()
        for line in f:
            line = line.strip('\n').split('\t')
            d[int(line[0])] = int(line[1])
    for key, value in d.items():
        dist[key].append(value)
with open('beadDist100iter.tsv', 'w') as f:
    for i in range(1,7,1):
        print(i)
        line = '\t'.join([str(x) for x in dist[i]])
        f.write('{0}\t{1}\n'.format(i, line))