#!/bin/bash

# Set sample name
# SAMPLE=GDIA_GAD_1_S4

SAMPLE=$(echo ${1} | sed "s/_R1_\001\.fastq\.gz//")

# Trim 5' sequencing primer
# -g ^ADAPTER 5' anchored sequence
# 2> redirects command line output to a txt file

cutadapt -g ^N{3}CCGGGGACTTATCATCCAACCTGT \
    ${SAMPLE}_R1_001.fastq.gz > \
    ../processed-seqs/${SAMPLE}_tn_trimmed.fastq.gz 2> ../processed-seqs/${SAMPLE}_trimming-info.txt \
    --discard-untrimmed

# change directory

cd ~/gdj-tnseq/processed-seqs/

# Use Bioawk to remove the 3' trailing sequence after MmeI cutsite
# I chose to do a short read with a perfect match. 
# One could choose to change this step to include a longer read and allow mismatches. 

~/sw/bioawk/bioawk -c fastx \
    '{ print "@"$name" "$comment; print substr($seq,1,16); print "+"; print substr($qual,1,16); }' \
    ${SAMPLE}_tn_trimmed.fastq.gz > ${SAMPLE}_16bp.fastq

# Separate seqs that begin with TA

~/sw/bioawk/bioawk -c fastx \
    '{ if ($seq ~ /^TA/) { print "@"$name" "$comment; print $seq; print "+"; print $qual; }}' \
    ${SAMPLE}_16bp.fastq > ../genome-info/${SAMPLE}_16bp_TAonly.fastq

# Move for alignment 

cd ~/gdj-tnseq/genome-info/

# Align TAonly seqs to reference genome
# -v 0 Report alignments with at most 0 mismatches
# -q Input files are .fastq (default)
# -m 1 Suppress all alignments if more than 1 reportable alignments exists
# -S Prints a header with @HD, @SQ and @PG lines to the .sam output files

~/sw/bowtie-1.2.3/bowtie -v 0 -q -m 1 -S gdj \
    ${SAMPLE}_16bp_TAonly.fastq \
    ../processed-seqs/${SAMPLE}_TAonly.sam 2> ../processed-seqs/${SAMPLE}_alignment-info.txt

# Move file when done with alignment

mv ${SAMPLE}_16bp_TAonly.fastq ~/gdj-tnseq/processed-seqs

cd ~/gdj-tnseq/processed-seqs/

# Count seqs aligning to TA sites

python3 ../scripts/count_TAs.py \
    ../genome-info/TA_sites_unfiltered.txt \
    ${SAMPLE}_TAonly.sam \
    > ../final-counts_tsv/${SAMPLE}_TAonly_counts.txt
