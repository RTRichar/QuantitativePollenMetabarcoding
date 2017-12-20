#!/usr/bin/env python

#this script takes and input taxonomy file, sys.argv[1], and revises Genbank lineages which have labels at the leaves of the tree (family, genus or species) but are 'undefined' at intermediate positions of the tree (class, order or family). To add labels to these entries, the lowest resolution labelled rank is used as a tag for the intermediate nodes (e.g 'o__Asterales, f__Undef, g__Solidago' becomes 'o__Asterales, f__gen_Solidago, g__Solidago'). This is not a great solution, however, given the nature of hierarchical sequence classification, it is necessary to do in order to provide the classifier with distinct lineages for each sequence. Otherwise, the classifier would not be able to destinguish a f__undef Solidago linneage from and f__undef Ranuculus lineage.

import sys
import re

ifile = open(sys.argv[1], 'r')

Raw_Taxonomies = {}
for line in ifile:
        Gi = line.split('\t')[0]
        Tax = line.split('\t')[1].strip('\n')[:-1]
        lin = Tax.split(';')
        Raw_Taxonomies[Gi] = lin

#print Raw_Taxonomies, '\n\n'

for i in range(8):
        for key in Raw_Taxonomies:
                lin = []
                for i in range(-1, -1 - len(Raw_Taxonomies[key]), -1):
                        if re.search( r'undef_', Raw_Taxonomies[key][i]) is not None:
                                lin.append(str('urs_' + Raw_Taxonomies[key][i + 1]))
                        else:
                                lin.append(str(Raw_Taxonomies[key][i]))
                Raw_Taxonomies[key] = lin[::-1]

#print Raw_Taxonomies

for key in Raw_Taxonomies:
        lin = Raw_Taxonomies[key]
        lintwo = []
        for i in range(len(Raw_Taxonomies[key])):
                if i == 0:
                        pre = 'k__'
                if i == 1:
                        pre = 'p__'
                if i == 2:
                        pre = 'c__'
                if i == 3:
                        pre = 'o__'
                if i == 4:
                        pre = 'f__'
                if i == 5:
                        pre = 'g__'
                if i == 6:
                        pre = 's__'
                if re.search( r'^[a-z]__', Raw_Taxonomies[key][i]) is None:
                        lintwo.append(str(pre + Raw_Taxonomies[key][i]))
                else:
                        lintwo.append(str(Raw_Taxonomies[key][i]))
        Raw_Taxonomies[key] = lintwo

for key, value in Raw_Taxonomies.iteritems():
        myString = str(";".join(value) + ';' + '\n')
        Raw_Taxonomies[key] = myString

### write final version of dictionary to file
outfile = open(str(str(sys.argv[1]).split('.')[0] + '_revised' + '.tax'), 'w')
for key, value in Raw_Taxonomies.items():
        outfile.write( str(key) + '\t' + str(value).strip('[ | ]'))

outfile.close()
ifile.close()
