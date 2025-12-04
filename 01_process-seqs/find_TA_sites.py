#!/usr/bin/env python3

## This script finds and reports the position of all TA sites in a genome
## input: python3 find_TA_sites.py argv[1] > output_file.txt
## argv[1] is the genome fasta sequence of the organism
## The output file has 1 column: 1-based index of the T in TA on the Fwd (+) strand

import sys

fasta = [] # use an array for effecient concatenation of the lines

with open(sys.argv[1]) as fp:
    for line in fp:
        # ignore the first line if it starts with >
        if line[0] == ">":
            first_line = line.split(' ')
            chrm = first_line[0][1:]
        else:
            # now only take the portion of the line up until a tab character
            split = line.split('\t')
            read = split[0].rstrip()
            fasta.append(read)

# join the lines together into a single string, and uppercase it
genome = "".join(fasta).upper()

index = genome.find('TA')
while (index != -1):
    print ("%s\t%d" %(chrm, index + 1))
    index = genome.find('TA', index + 2)
