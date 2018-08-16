#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 11:25:52 2018

@author: Zewei Song
@email: songzewei@genomics.cn
"""
from __future__ import print_function
from __future__ import division
import argparse
from metaSeq import io as seqIO
import json
import textwrap

parser = argparse.ArgumentParser()
parser.add_argument('-i', help='Input beadJson format file.')
parser.add_argument('-o', help='The output FASTQ file. Barcode will be put at -2, while /1 and /2 at -1.')
parser.add_argument('-bp', default=-1, help='The position of barcode in beadJson. usually -1 or -2')
args = parser.parse_args()