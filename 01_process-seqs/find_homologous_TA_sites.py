#!/usr/bin/env python3

## This script finds and reports replicate sequences within the file
## input example: python3 find_homologous_TA_sites.py argv[1] > output_file.txt
## argv[1] is the genome fasta file of the organism
## The output file has 2 columns: 1) the non-unique TA sequence; 2) 1-based index of the T in TA on the Fwd (+) strand

import sys

homology = {}

with open(sys.argv[1]) as fp:
	for line in fp:
		split = line.split('\t')
		sequence = split[0]
		list = homology.get(sequence)
		if list == None:
			list = []
			homology[sequence] = list
		list.append(line)

for sequence, list in homology.items():
	if len(list) > 1:
		for line in list:
			print (line.rstrip())
