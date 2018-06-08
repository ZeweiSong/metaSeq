#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  8 15:40:00 2018

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%%
from __future__ import print_function
from __future__ import division
from metaSeq import io as seqIO
from metaSeq import qc as seqQC
import json
import argparse
import textwrap
import os
import time

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
                                        Do QC per bead. Discard reads with expected error >= 1.
                                        For paired reads, both need to pass QC to stay.
                                        Output to a new JSON bead file.
                                        ------------------------
                                        By Zewei Song
                                        Environmental Ecology Lab
                                        Institute of Metagenomics
                                        BGI-Research
                                        songzewei@genomics.cn
                                        songzewei@outlook.com
                                        ------------------------'''))
parser.add_argument('-i', help='Input JSON bead file.')
parser.add_argument('-o', help='Output JSON bead file.')
args = parser.parse_args()

beadFile = args.i
outputFile = args.o
logFile = 'stlft_bead_qc.log'
if os.path.isfile(logFile):
    os.remove(logFile)
fl = open(logFile, 'a')

p_dict = seqQC.qual_score()
t1 = time.time()
with open(outputFile, 'w') as f:
    count = 0
    keep = 0
    for bead in seqIO.bead_json(beadFile):
        count += 1
        if count % 1000000 == 0:
            t2 = time.time()
            fl.write('Processed {0} beads, kept {1} ({2:3.2f}%), used {3:3.1f}s.'.format(count, keep, keep/count*100, t2-t1))
        filterRecords = []
        for record in list(bead.values())[0]:
            if len(record) == 2:
                if seqQC.ee(record[1], p_dict) < 1:
                    filterRecords.append([record[0]])
            else:
                if seqQC.ee(record[1], p_dict) < 1 and seqQC.ee(record[3], p_dict) < 1:
                    filterRecords.append([record[0], record[2]])
        if len(filterRecords) > 0:
            keep += 1
            f.write('%s\n' % json.dumps({list(bead.keys())[0]:filterRecords}))
t2 = time.time()
fl.write('Processed {0} beads, kept {1} ({2:3.2f}%), used {3:3.1f}s.'.format(count, keep, keep/count*100, t2-t1))
fl.close()