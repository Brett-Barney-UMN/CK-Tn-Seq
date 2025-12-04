## This script assigns gene names to TA_sites.txt
## argv[1] is a list of TA sites
## arvg[2] is a gene list
## The output file has two columns: 1) Index; 2) Gene

import sys

# Add TA_sites to dictionary

TA_dict = {}

with open(sys.argv[1]) as fp:
    for line in fp:
        split = line.split() # splits into columns--treates consecutive spaces as 1 single tab
        if split[0] != '': # if line is not empty
            TA_site = int(split[1]) # if yes, get TA index
            TA_dict[TA_site] = [0] # add TA index to dict, initialize with null value


# Add genes to dictionary

gene_dict = {}; prev_end = 0

for line in open(sys.argv[2]):
    split_line = line.split('\t')
    chrm = split_line[0]; IG = 'IG_' + chrm
    gene = split_line[1]; start = int(split_line[4]); end = int(split_line[5])
    if gene in gene_dict:
        old_start = int(gene_dict[gene][0])
        old_end = int(gene_dict[gene][1])
        if start < old_start:
            gene_dict[gene][0] = start
        if end > old_end:
            gene_dict[gene][1] = end
    else:
        gene_dict[gene] = [start, end, chrm]
    IG = 'IG_' + gene
    if start > prev_end + 1:
        gene_dict[IG] = [prev_end + 1, start - 1, chrm]
    prev_end = end

gene_dict['IG_' + chrm + 'end'] = [prev_end + 1, int(sys.argv[3]), chrm] # last gene nt, last nt, chrm

# Add genes to TA_dict

for k, v in gene_dict.items():
    start = v[0]; end = v[1]; chrm = v[2]
    for i in range(start, end + 1): # add + 1 to include the end point
        if i in TA_dict:
            TA_dict[i].append(k)

# Print

TA_sites = list(TA_dict.keys())
TA_sites.sort()

for i in TA_sites:
    print (chrm, i, TA_dict[i][1], sep = "\t")
