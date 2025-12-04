#!/usr/bin/env python3

## This script parses a gff3 file into a simpler version
## Input: python3 parse_gff_file.py argv[1] > output_file.txt
## argv[1] is the gff file to parse
## Output file has 6 columns:
## 1) chr; 2) locus_tag; 3) name; 4) description; 5) start; 6) end

import sys

with open(sys.argv[1]) as fp:
    first_line = fp.readline()
    accession_num = first_line.split()[1]
    for line in fp:
        if line.strip() != "": # ignore empty lines
            if line[0][0] != "#": # ignore lines starting with #
                split = line.split(None, 8)
                if split[2] == "CDS":

                    start = split[3]
                    end = split[4]
                    strand = split[6]
                    attributes = split[8].split(';')

                    locus_tag = [i for i in attributes if i.startswith('locus_tag')]
                    product = [i for i in attributes if i.startswith('product')]
                    gene = [i for i in attributes if i.startswith('gene')]

                    locus_tag = locus_tag[0][len("locus_tag="):]
                    product = product[0][len("product="):]

                    if not gene:
                        gene = locus_tag
                    else:
                        gene = gene[0][len("gene="):]

                    print (accession_num, locus_tag, gene, product, start, end, strand, sep = '\t')
