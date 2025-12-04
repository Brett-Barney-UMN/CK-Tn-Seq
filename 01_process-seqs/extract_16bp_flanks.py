#!/usr/bin/env python3

## This script extracts the expected 16bp flanks associated with each TA site
## input: python3 extract_16bp_flanks.py argv[1] > output_file.txt
## argv[1] is the genome fasta file of the organism
## The output file has 2 columns:
## 1) TA + 14 nucleotides downstream of TA;
## 2) 1-based index of the T in TA on the Fwd (+) strand

import sys

fasta = [] # use an array for effecient concatenation of the lines

with open(sys.argv[1]) as fp:
    for line in fp:
        # ignore the first line if it starts with >
        if line[0] != ">":
            # now only take the portion of the line up until a tab character
            split = line.split('\t')
            read = split[0].rstrip()
            fasta.append(read)

# join the lines together into a single string, and uppercase it
genome = "".join(fasta).upper()

index = genome.find('TA')
while (index != -1):
    downstream_sequence = genome[index:index + 16]
    print ("{}\t{}".format(downstream_sequence, (index + 1)))
    index = genome.find('TA', index + 2)
