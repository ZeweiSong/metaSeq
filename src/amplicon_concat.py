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
import json
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', type=argparse.FileType('r'), nargs='+', default=sys.stdin, help='profile list separated by space')
parser.add_argument('-biom_out', default='combined.biom', help='Output biom file name')
args = parser.parse_args()

biomFile = args.biom_out

biomList = []
for f in args.input:
    biomList.append(f)
if len(biomList) <= 1:
    print('Found only one biom profile, will still give you a (brand) new one.')
    biomProfile = Table.from_json(json.load(biomList[0]))
    with open(biomFile, 'w') as f:
        biomProfile.to_json('Generated_by_almighty_metaSeq', f)
else:
    print('Found {0} biom profiles under.'.format(len(biomList)))
    biomProfile = Table.from_json(json.load(biomList[0]))

    for f in biomList[1:]:
        biomProfile = biomProfile.concat([Table.from_json(json.load(f))])

    with open(biomFile, 'w') as f:
        biomProfile = biomProfile.sort()
        biomProfile.to_json('Generated_by_almighty_metaSeq', f)
    print('Concatenated {0} profiles into {1}.'.format(len(biomList), biomFile))