#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 11:23:14 2018

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%%
from __future__ import print_function
from __future__ import division
#import biom.table as Table
import biom

# Load table
t = biom.load_table('combined.taxa.biom')

#%%
print(t) # Print table
print(t.ids()) # Print sample names
print(t.ids(axis='observation')) # Print OTU names
print(t.nnz) # Number of non zero entries
print(t.head()) # Print head

#%%
t = t.sort() # Sort by sample names
print(t.sum(axis='sample'))
print(t.sum(axis='observation'))

#%%
# Sort the table by the ave/sum abundance of observations/OTUs
ids = [x for _,x in sorted(zip(t.sum(axis='observation'),t.ids(axis='observation')), reverse=True)]
# print(ids)
t_sorted = t.sort_order(ids, axis='observation')
print(t_sorted)

#%%
# Locate the abundace of a sample and an OTU
t[t.index('Fusarium_oxysporum|SH200825.07FU', 'observation'),t.index('502', 'sample')]
# OR
t.get_value_by_ids('Fusarium_oxysporum|SH200825.07FU', '502')