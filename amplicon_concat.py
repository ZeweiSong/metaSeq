#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 10:36:47 2018

Concat all biom OTU table into a mighty one.

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%%
from __future__ import print_function
from __future__ import division
import argparse
from biom.table import Table
import os

parser = argparse.ArgumentParser()
parser.add_argument('-if', '--inputFolder', help='path/to/the/profile/folder')
parser.add_argument('-biom_out', default='combined.biom', help='Output biom file name')
#parser.add_argument('-table_out', default='combined.profile', help='Outpuit biom file name')
args = parser.parse_args()

inputFolder = args.inputFolder
biomFile = args.biom_out
#table_out = args.table_out

fileList = os.listdir(inputFolder)
biomList = [i for i in fileList if i[-5:] == 'biom']