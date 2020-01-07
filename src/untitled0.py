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
from math import lgamma
from math import exp
#print(ceil(lgamma(100000+1)/log(10)))
#print(lgamma(100000+1))

def occur_prob(n, m, k):
    return 1 - exp(lgamma(n-m) - lgamma(n) + lgamma(n-k) - lgamma(n-k-m))

n = 10**5
m = n * 0.001
k = 10**4

print(occur_prob(n, m, k))

#%% Draw the curve
n = 10**5
k = 10**4
p = [i/n for i in range(1, n//20, 10)]

line = []
for item in p:
    xy = (item*100, occur_prob(n, int(n*item), k))
    line.append(xy)

#%%
n = 10**8//1
k = [50000,10000, 5000, 2000, 1000]
p = [i/n for i in range(1, n//100, 1)]

line = {}
for value_k in k:
    line[value_k] = [(i*1000, occur_prob(n, int(n*i), value_k)) for i in p]

import matplotlib.pyplot as plt
for key, value in line.items():
    print(key)
    plt.plot([i[0] for i in value], [i[1] for i in value])
#%%
import random
import bisect
# Read in a biom otu table
# Convert to relative abundance version 

# n total sum
# m count of feature i
# k sampling depth
def occur_prob(n, m, k):
    return 1 - exp(lgamma(n-m) - lgamma(n) + lgamma(n-k) - lgamma(n-k-m))

# SET
k = 10000 # Sampling depth
p_thre = 0.97

# Read in biom table
# pick one sample as dictionary
sample = {}
n = sum([i for i in sample.values()]) # The sum of the sample
p = [(i, occur_prob(n, j, k)) for i, j in sample.items()] # Store the key as a list of lists.
p.sort(key=lambda i:i[1], reverse=True)

ra = [(i[0], sample[i[0]]/n*100) for i in p if i[1] > p_thre] # Relative abundance for species not entering the drawing
used_p = sum([i[1] for i in ra]) # Used proportion of k

# Get the left abundance from k
k_leftover = k * (1 - used_p)
sample_tail = [(i[0], sample[i[0]]) for i in p if i[1] <= p_thre]
# Randomly draw i times, pick the one with median length
itertime = 1000

running_total = []
for index, item in enumerate(sample_tail):
    running_total += [index] * item[1]

sample_pool = []
for i in range(itertime):
    working_sample = random.sample(running_total, k_leftover)
    names = tuple(set(working_sample))
    left = {i:0 for i in names}
    for item in working_sample:
        left[i] = left.get(i, 0) + 1
    sample_pool.append({sample_tail[i]:j for i,j in left.items()})