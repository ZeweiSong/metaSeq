#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 10:19:30 2018

@author: Zewei Song
@email: songzewei@genomics.cn
"""
from __future__ import print_function
from __future__ import division
import argparse
from metaSeq import amplicon
from biom.table import Table
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-i', nargs='+', help='File names of the two alignment output.')
parser.add_argument('-t', '--tab_out', help='File name of tab-delimited format output.')
parser.add_argument('-b', '--biom_out', help='File name of biom format output.')
method = parser.add_mutually_exclusive_group(required=True)
method.add_argument('-g', '--greedy', action='store_true', default=True, help='Option to enable greedy method. [Default is greedy capitalists]')
method.add_argument('-l', '--less_greedy', action='store_true', default=False, help='Option to enable less greedy method, equally greedy is nTargets = 1. [Enable me to disable the greedy capitalists]')
parser.add_argument('-sn', default='SampleX', help="Sample name if specified.")
score = parser.add_mutually_exclusive_group()
#s = score.add_argument_group()
score.add_argument('-ec', action='store_true', default=True, help='Specify to score using Effective Count [Default]')
score.add_argument('-median', action='store_true', default=False, help='Specify to score using median of abundance.')
score.add_argument('-ave', action='store_true', default=False, help='Specify to score using average of abundance.')
args = parser.parse_args()

alnString = []
for item in args.i:
    alnString.append(item)
targetNumber = len(alnString) # The number of targets used (normally this is 2, e.g. ITS1 forward & ITS1 reverse)
sampleName = args.sn
tabFile = args.tab_out
biomFile = args.biom_out

greedy = args.greedy
less_greedy = args.less_greedy
if less_greedy:
    greedy = False

ec = args.ec
median = args.median
ave = args.ave
if median:
    weight = 'median'
elif ave:
    weight = 'average'
elif ec:
    weight = 'ec'

print('Ah the bloody capitalism came to {0}!'.format(sampleName))
mode = 'greedy'
if not greedy:
    mode = 'less greedy'
print('Released are the {0} capitalists.'.format(mode))
if not greedy:
    print("\tI'll spit out some profit to miximize mine.\n")
else:
    print('\tAll profit is mine!\n')
alnNormalized = amplicon.initAlignment(alnString)
wtaProfile = amplicon.winnerTakeAll(alnNormalized, progress=True, greedy=greedy, weight=weight)
print('Winner take all profile found! {0} references survived.'.format(len(wtaProfile)))
profile = list(wtaProfile.items())
profile.sort(key=lambda x:x[1], reverse=True)

print('Tab-delimited profile wrote to {0}.'.format(tabFile))
with open(tabFile, 'w') as f:
    f.write('{0}\t{1}\n'.format('Reference', sampleName))
    for item in profile:
        if int(item[1]) >= 1:
            f.write('{0}\t{1}\n'.format(item[0], int(item[1])))

# Output as biom format
biomTable = Table(np.array([[i[1]] for i in profile]), [i[0] for i in profile], [sampleName], table_id='single_sample_biom')
print('Biom JSON format profile wrote to {0}.'.format(biomFile))
with open(biomFile, 'w') as f:
    biomTable.to_json('Generated_by_metaSeq', f)