#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 16:59:32 2017

Test on download file in Python

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%%
from __future__ import print_function
from __future__ import division
from urllib import request

request.urlretrieve('ftp://ftp.ncbi.nih.gov/genomes/refseq/fungi/assembly_summary.txt', 'assembly_summary_test.txt')
