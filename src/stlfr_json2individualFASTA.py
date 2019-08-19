#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 09:59:28 2018

@author: Zewei Song
@email: songzewei@genomics.cn
"""
from __future__ import print_function
from __future__ import division
from metaSeq import io as seqIO
from metaSeq import bead
from metaSeq import kmer
import json
from itertools import combinations
import textwrap
import argparse
import re

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
                                        Save sequence individually into each bead. Open the script to see a detail document.
                                        ------------------------
                                        By Zewei Song
                                        Environmental Ecology Lab
                                        Institute of Metagenomics
                                        BGI-Research
                                        songzewei@genomics.cn
                                        songzewei@outlook.com
                                        ------------------------'''))
parser.add_argument('-i', help='merged fasta file.')
parser.add_argument('-d', help='Output directory.')
parser.add_argument('-z', help='save fa with missing tag(0000) into a seperated directory.')

args = parser.parse_args()
#%% Try to sample bead with more than 1 fragments
inputBeadJson = args.i

beadPool = []
with open(inputBeadJson, 'r') as f:
    for line in f:
        beadPool.append(bead.beadSequence(json.loads(line)))
print('Find {0} beads.'.format(len(beadPool)))

#%% Get beads with more than 1 fragment
beadPool_1 = []
for item in beadPool:
    if len(item.fragments) > 1:
        beadPool_1.append(item)
print('Find {0} bead with more than 1 fragments.'.format(len(beadPool_1)))

#%% Write these bead to FASTA files
outputFolder = args.d
code0Folder = args.z

for item in beadPool_1:
    if re.match("0000",item.barcode,0) :
        item.fastaWrite(folder=code0Folder)
    else:
        item.fastaWrite(folder=outputFolder)
