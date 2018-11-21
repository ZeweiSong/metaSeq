#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 11:42:28 2018

A place holder for split function.

@author: Zewei Song
@email: songzewei@genomics.cn
"""
from __future__ import print_function
from __future__ import division
import argparse
from metaSeq import io as seqIO

parser = argparse.ArgumentParser()
parser.add_argument('-r1', help='R1 reads')
parser.add_argument('-r2', help='R2 reads if needed')
parser.add_argument('-b', '--barcode', help='Barcode mapping file with sample info')
parser.add_argument('-o', '--output', help='Output folder for splitted files')
args = parser.parse_args()
r1File = args.r1
r1f = seqIO.sequence(r1File, fastx='q', gz=True)
if args.r2:
    r2File = args.r2