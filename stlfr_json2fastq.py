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

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
                                        Convert the pair end merged FASTQ files into JSON format.
                                        The JSON file can be read in using the Class bead_json().
                                        By Zewei Song
                                        Environmental Ecology Lab
                                        Institute of Metagenomics
                                        BGI-Research
                                        songzewei@genomics.cn
                                        songzewei@outlook.com
                                        ------------------------'''))
parser.add_argument('-a', help='FASTQ file of assembled reads.')
parser.add_argument('-f', help='Unassembled forward reads.')
parser.add_argument('-r', help='Unassembled reverse reads.')
parser.add_argument('-o', help='Output bead file in JSON format.')
args = parser.parse_args()