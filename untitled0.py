#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 15:46:46 2018

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%%
from __future__ import print_function
from __future__ import division
sra = []
with open('SRR_Acc_List-1.txt', 'r') as f:
    for line in f:
        sra.append(line.strip('\n'))
with open('sra_down.sh', 'w') as f:
    for line in sra:
        f.write('wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/{0}/{1}/{1}.sra'.format(line[:3], line))