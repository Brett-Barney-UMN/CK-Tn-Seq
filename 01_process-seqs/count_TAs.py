# This script takes each read in the SAM file and tallies them to each TA site in the genome.

# input:  python3 count_TAs.py argv[1] argv[2] > output.txt
## argv[1] list of annotated TA sites (unfiltered)
## arvg[2] mapped reads in SAM format (single end mapping)

# The output file has 6 columns:
## 1) name of chromosome;
## 2) start position of TA site in genome;
## 3) the locus_tag in which each TA site sits;
## 4) reads from insertions at that TA site on the Forward (+) strand;
## 5) reads from insertions at that TA site on the Reverse (-) strand

import sys

TA_counts = {} # create dictionary, will be nested in future steps

with open(sys.argv[1]) as TA_sites:
    for line in TA_sites: # sys.argv[1]
        split = line.split() # splits string into columns--treates consecutive spaces as 1 single tab
        if line.strip() != "": # ignore empty lines
            chrm = str(split[0]) # get chrm
            TA_site = int(split[1]) # get TA index
            locus_tag = split[2]
            if chrm in TA_counts:
                TA_counts[chrm].update({TA_site: [0, 0, locus_tag]})
            else:
                TA_counts[chrm] = {TA_site: [0, 0, locus_tag]}

# count reads at each TA site from SAM output

with open(sys.argv[2]) as map:
    for line in map:  # sys.argv[2]
        split = line.split('\t')
        if split[0][0] != "@": # ignores headers (analyzes rows don't have an @ symbol)
            chrm = str(split[2])
            if split[1] == '0':  # if mapping is to the plus strand
                TA_site = int(split[3])
                TA_counts[chrm][TA_site][0] += 1
            if split[1] == '16': # if mapping is on the minus strand
                TA_site = int(split[3]) - 2 + len(split[9])
                TA_counts[chrm][TA_site][1] += 1 # if read has been found before, add 1

# print to file

for chrms, site_info in TA_counts.items():
    for i in site_info:
        print(chrms, i, TA_counts[chrms][i][2], TA_counts[chrms][i][0], TA_counts[chrms][i][1], sep = "\t")
