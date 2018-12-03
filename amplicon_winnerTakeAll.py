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
parser.add_argument('-sn', default='SampleX', help="Sample name if specified.")
parser.add_argument('-suppress_progress', action='store_true', help='Indicators for suppressing progress.')
args = parser.parse_args()

alnString = []
for item in args.i:
    alnString.append(item)
targetNumber = len(alnString) # The number of targets used (normally this is 2, e.g. ITS1 forward & ITS1 reverse)
sampleName = args.sn
tabFile = args.tab_out
biomFile = args.biom_out
progress = True
if args.suppress_progress:
    progress = False

print('Found {0} targets from your input.'.format(targetNumber))
print('Normalizing alignments ...')
alnNormalized = amplicon.initAlignment(alnString)
print('Only one winner survives each round ...')
wtaProfile = amplicon.winnerTakeAll(alnNormalized, progress=progress)
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