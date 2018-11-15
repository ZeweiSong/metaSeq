#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 16:58:07 2018

Convert the GreenGene taxonomy into Burst format

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%%
from __future__ import print_function
from __future__ import division
inputFile = '61_otu_taxonomy.txt'
outputFile = inputFile.replace('_taxonomy.txt', '.tax')