#!/usr/bin/env python3

## This script filters out nonpermissive and homologous sites

import sys

# Add TA_sites to arrays

TA_sites = {}

with open(sys.argv[1]) as fp:
    for line in fp:
        split_line = line.split()
        chrm = split_line[0]
        TA_site = int(split_line[1])
        locus_tag = split_line[2]
        TA_sites[TA_site] = [chrm, locus_tag]

nonpermissive_sites = []

with open(sys.argv[2]) as fp:
    for line in fp:
        split_line = line.split()
        TA_site = int(split_line[1])
        nonpermissive_sites.append(TA_site)

homo_sites = []

with open(sys.argv[3]) as fp:
    for line in fp:
        split_line = line.split()
        TA_site = int(split_line[1])
        homo_sites.append(TA_site)

for i in nonpermissive_sites:
     if i in TA_sites:
         del TA_sites[i]

for i in homo_sites:
     if i in TA_sites:
          del TA_sites[i]

# print
    
for i in TA_sites:
    print(TA_sites[i][0], i, TA_sites[i][1], sep = "\t")
