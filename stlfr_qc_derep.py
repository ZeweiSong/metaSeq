#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 15:42:25 2018

@author: Zewei Song
@email: songzewei@genomics.cn
"""

from __future__ import print_function
from __future__ import division
import json
from metaSeq import io as seqIO
from metaSeq import qc as seqQC
from metaSeq import bead

inputFile = 'mock.json'
maxee = 1
outputFile = 'mock.qc.derep.json'

mockRaw = []
with open(inputFile, 'r') as f:
    for line in f:
        mockRaw.append(json.loads(line))