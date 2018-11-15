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

parser = argparse.ArgumentParser()
parser.add_argument('-i', help='File names of the two alignment output.')
parser.add_argument('-t', '--tab_out', help='File name of tab-delimited format output.')
#parser.add_argument('-b', '--biom_out', help='File name of biom format output.')
parser.add_argument('-sn', default= 'Sample', help="Sample name if specified.")
args = parser.parse_args()

alnString = args.i
#alnString = 'even1.forward.unite_ee.burstFORAGE97.b6,even1.reverse.unite_ee.burstFORAGE97.b6' # Provide alignments files on different targets
alnString = alnString.strip('\n').split(',')
TargetNumber = len(alnString) # The number of targets used (normally this is 2, e.g. ITS1 forward & ITS1 reverse)
sampleName = args.sn
tabFile = args.tab_out

alnNormalized = amplicon.initAlignment(alnString)
wtaProfile = amplicon.winnerTakeAll(alnNormalized)
profile = list(wtaProfile.keys())
profile.sort(key=lambda x:x[1], reverse=True)
with open(tabFile, 'w') as f:
    f.write('{0}\t{1}\n'.format('Reference', sampleName))
    for item in profile:
        if item[1] >= 1:
            f.write('{0}\t{1}\n'.format(item[0], int(item[1])))

# TODO need to add output report
# Add some output