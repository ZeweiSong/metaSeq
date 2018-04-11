#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 16:59:32 2017

@author: songz
"""
from __future__ import print_function
from __future__ import division

#%% Generate a look up dictionary for quality score, the value is P instead of Q
def qual_score():
#    qual = {}
    qual_list = ['"', '#', '$', '%', '&', "'", '(', ')', '*', '+', ',', '-', '.', '/', \
                '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', ':', ';', '<', '=', \
                '>', '?', '@', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']
#    i = 1
#    for char in qual_list:
#        qual[i] = char
#        i += 1
    p = {}
    i = 1
    for char in qual_list:
        p[char] = 10**(-1*i/10)
        i += 1
    p['!'] = 1 # Add ! (means N)
    return p


#%% Trim a FASTQ record from the end, until MaxError/base lower than the threshold
def trunc_ee_rate(record, p_dict, rate = 0.01):
    # Convert Q score to P
    p = []
    for qchar in record[3]:
        p.append(p_dict[qchar])
    
    er_sum = sum(p)
    len_record = len(record[1])
    trim_pos = len_record
#    print(er_sum)
#    print(er_sum/trim_pos)
#    print(rate)
#    print(trim_pos)
#    print(er_sum/trim_pos)
#    print(er_sum/trim_pos >= rate)
    # Check from the last base
    while er_sum/trim_pos >= rate and trim_pos > 1:
        trim_pos -= 1
        er_sum = sum(p[0:trim_pos])
        #print(er_sum, trim_pos)        
    return [record[0], record[1][0:trim_pos], record[2], record[3][0:trim_pos]]
    #for i in range(len_record, 0, -1):