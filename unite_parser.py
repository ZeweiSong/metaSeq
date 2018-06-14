#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 16:59:32 2017

@author: songz
"""
#%% Read in UNITE
from __future__ import print_function
from __future__ import division
from metaSeq import io
unite_file = 'sh_general_release_dynamic_01.12.2017.fasta'
unite = []
for record in io.sequence(unite_file):
    unite.append(record)
levels = {'k':[], 'p':[], 'c':[], 'o':[], 'f':[], 'g':[], 's':[], 'sh':[]}
level = ['k', 'p', 'c', 'o', 'f', 'g', 's', 'sh']
#%% Split UNITE taxonomy
unite_taxon = []
for record in unite:
    taxon = record[0].split('|')
    taxon = taxon[0:4] + taxon[4].split(';') + [taxon[2]]
    unite_taxon.append(taxon)

for record in unite_taxon:
    for i in range(4,12,1):
        if record[i] not in levels[level[i-4]]:
            levels[level[i-4]].append(record[i])
#%% Phylum count under each Kingdom
phylum_count = {}
for key in levels['k']:
    phylum_count[key] = []
for record in unite_taxon:
    if record[5] not in phylum_count[record[4]]:
        phylum_count[record[4]].append(record[5])
class_count = {}
for key in levels['p']:
    class_count[key] = []
for record in unite_taxon:
    if record[6] not in class_count[record[5]]:
        class_count[record[5]].append(record[6])
#%%
unite_dict = {}
for record in unite_taxon:
    unite_dict[record[4]] = {}
for record in unite_taxon:
    unite_dict[record[4]][record[5]] = {}
for record in unite_taxon:
    unite_dict[record[4]][record[5]][record[6]] = {}
for record in unite_taxon:
    unite_dict[record[4]][record[5]][record[6]][record[7]] = {}
for record in unite_taxon:
    unite_dict[record[4]][record[5]][record[6]][record[7]][record[8]] = {}
for record in unite_taxon:
    unite_dict[record[4]][record[5]][record[6]][record[7]][record[8]][record[9]] = {}
for record in unite_taxon:
    unite_dict[record[4]][record[5]][record[6]][record[7]][record[8]][record[9]][record[10]] = {}
for record in unite_taxon:
    unite_dict[record[4]][record[5]][record[6]][record[7]][record[8]][record[9]][record[10]][record[11]] = {}
#%%
print('There are {0} Kingdoms.'.format(len(unite_dict)))
for k, t in unite_dict.items():
    print('{0} has {1} Phylums.'.format(k, len(t)))
    for p, t2 in t.items():
        print('\t{0} has {1} Class.'.format(p, len(t2)))
        for c, t3 in t2.items():
            print('\t\t{0} has {1} Order.'.format(c, len(t3)))
            for o, t4 in t3.items():
                print('\t\t\t{0} has {1} Family.'.format(o, len(t4)))
                for f, t5 in t4.items():
                    print('\t\t\t\t{0} has {1} Genus.'.format(f, len(t5)))
                    for g, t6 in t5.items():
                        print('\t\t\t\t\t{0} has {1} Specis.'.format(g, len(t6)))
#%% Write to XML format
fc = open('report_count.txt', 'w')
fc.write('<Life attribute={0} child={1}>\n'.format('Life', len(unite_dict)))
#fc.write('There are {0} Kingdoms.\n'.format(len(unite_dict)))
for k, t in unite_dict.items():
    fc.write('<kingdom attribute={0} child={1}>\n'.format(k, len(t)))
    for p, t2 in t.items():
        fc.write('\t<phylum attribute={0} child={1}>\n'.format(p, len(t2)))
        for c, t3 in t2.items():
            fc.write('\t\t<class attribute={0} child={1}>\n'.format(c, len(t3)))
            for o, t4 in t3.items():
                fc.write('\t\t\t<order attribute={0} child={1}>\n'.format(o, len(t4)))
                for f, t5 in t4.items():
                    fc.write('\t\t\t\t<family attribute={0} child={1}>\n'.format(f, len(t5)))
                    for g, t6 in t5.items():
                        fc.write('\t\t\t\t\t<genus attribute={0} child={1}>\n'.format(g, len(t6)))
                        for s, t7 in t6.items():
                            fc.write('\t'*6 + '<species attribute={0} child={1}>\n'.format(s, len(t7)))
                            for sh, t8 in t7.items():
                                pass
                                #fc.write('\t'*7 + '<SH attribute={0} child={1}>\n'.format(sh, len(t8)))
                                #fc.write('\t'*7 + '</SH>')
                            fc.write('\t'*6 + '</species>')
                        fc.write('\t\t\t\t\t</genus>\n')
                    fc.write('\t\t\t\t</family>\n')
                fc.write('\t\t\t</order>\n')
            fc.write('\t\t</class>\n')
        fc.write('\t</phylum>\n')
    fc.write('</kingdom>\n')
fc.write('</Life>\n')
fc.close()
#%% Generate key list for each level
levels = {'k':[], 'p':[], 'c':[], 'o':[], 'f':[], 'g':[], 's':[], 'sh':[]}
level = ['k', 'p', 'c', 'o', 'f', 'g', 's', 'sh']
for i in range(4, 11, 1):
    print(i)
    for record in unite_taxon:
        unite_dict[record[i]][record[i+1]] = {}