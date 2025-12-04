#!/bin/bash

# This script analyzes the TA sites for a provided chrm
# Outputs four files:
# 1 - list of all TA sites 
# 2 - list of nonpermissive TA sites 
# 3 - list of TA flanking seqs 
# 4 - Annotated list of unfiltered TA sites 

# 1. Set variables
fasta=${1}.fasta
rev_comp=${1}_revcomp.fasta
gff3=${1}.gff3

# 2. Find all TA sites
python3 find_TA_sites.py $fasta > ${1}_TA-sites.txt

# 3. Find just nonpermissive TA sites
python3 find_nonpermissive_TA_sites.py $fasta | sed "s/$/ \t${1}/" > ${1}_TA-sites_nonpermissive.txt

# 4. Find just homologous TA sites

# Extract 16bp flanks (TA[14bps])
python3 extract_16bp_flanks.py $fasta > 16bp_flanks_fwd.txt
python3 extract_16bp_flanks.py $rev_comp > 16bp_flanks_rev.txt

# Correct revcomp TA indices

# Write new file, last line first
tac 16bp_flanks_rev.txt > tmp.txt
# Extract the first column of sequences
cut -f1 tmp.txt > tmp2.txt
# Extract the second colum of indexes and combine with the sequences
cut -f2 16bp_flanks_fwd.txt | paste tmp2.txt - > tmp3.txt
# Combine fwd and rev sequences
cat tmp3.txt 16bp_flanks_fwd.txt > tmp4.txt
# Add chrm as column 
sed "s/$/ \t${1}/" tmp4.txt > ${1}_TA_flanks.txt

# Remove temp files
rm tmp*.txt
rm 16bp*.txt

# 5. Assign genes to TA sites

# get chrm length
chrm_length="$(head -1 $gff3 | cut -d" " -f4)"

# annotate TA sites
python3 annotate_TA_sites.py \
  ${1}_TA-sites.txt \
  ${1}_gene-list.txt \
  $chrm_length \
  > ${1}_TA-sites_unfiltered.txt
  
# 6. Remove unused files
rm ${1}_TA-sites.txt
  
