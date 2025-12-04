#!/bin/bash -l
#SBATCH --ntasks=2
#SBATCH --time=2:30:00
#SBATCH --mem=10g

cd ~/gdj-tnseq/raw-data_fastq

module load cutadapt

for i in *.fastq.gz;
do
  ~/gdj-tnseq/scripts/get-counts.sh ${i}
done
