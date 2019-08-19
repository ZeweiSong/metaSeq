#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 15:10:37 2018

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%%
from __future__ import print_function
from __future__ import division
import textwrap
import argparse
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
                                        Filter MASH file.
                                        ------------------------
                                        By Zewei Song
                                        Environmental Ecology Lab
                                        Institute of Metagenomics
                                        BGI-Research
                                        songzewei@genomics.cn
                                        songzewei@outlook.com
                                        ------------------------'''))
parser.add_argument('-i', help='the mash file.')
parser.add_argument('-o', help='the output dist file.')
parser.add_argument('-c', default=0.01, help='the threshold.')
args = parser.parse_args()

inputFile = args.i
threshold = float(args.c)
outputFile= args.o
pair = []
with open(inputFile, 'r') as f:
    for line in f:
        line = line.strip('\n').split('\t')
        source = line[0].split('/')[2]
        target = line[1].split('/')[2]
        distance = float(line[2])

        if source != target and 0.02 <= distance <= threshold:
            pair.append([source, target, distance])
with open(outputFile, 'w') as f:
    for line in pair:
        f.write('{0}\n'.format('\t'.join([str(i) for i in line])))
