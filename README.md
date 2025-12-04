# Tn-seq analysis 

There are three parts to this analysis: 

- [Part 1](https://github.umn.edu/umn-barney-lab/tn-seq/blob/master/README.md#part-1): Get TA sites (personal computer)
- [Part 2](https://github.umn.edu/umn-barney-lab/tn-seq/blob/master/README.md#part-2): Process sequences (MSI supercomputer)
- [Part 3](https://github.umn.edu/umn-barney-lab/tn-seq/blob/master/README.md#part-3): Calculate fitness (personal computer)
- [Part 4](https://github.umn.edu/umn-barney-lab/tn-seq/blob/master/README.md#part-4): Compare proteins using desktop Blast (personal computer)

Before getting started, you will need the following external data.
1. FASTA file (ex. CP001190.fasta)
2. Reverse complement of FASTA file (ex. CP001190_revcomp.fasta)
3. GFF3 file (ex. CP001190.gff3)

> Note: For the scripts to work, please follow the naming convention shown. Gdj data is being used as an example.

## How to download .fasta and .gff3 files

### Browser

1. Go to the KEGG website for your organism. For gdj > https://www.genome.jp/kegg-bin/show_organism?org=T00798.
1. Under Genome Information > Click Sequence. This link will bring you to the NCBI entry. For gdj you will have two links, CP001189 and CP001190.
3. On the NCBI website > Click Send to: Chose Complete File, File, and chose file format FASTA. *file #1*
4. Then > Customize View > Show reverse complement > update view and repeat the steps above. *file #2*
5. On the NCBI website > Click Send to: Chose Complete File, File, and chose file format GFF3. *file #3*
6. Repeat steps for each chromosome and/or plasmid in your organism.
7. Place all files in your project folder and record the date and how each file was downloaded. Ex. Downloaded via web browser from https://www.ncbi.nlm.nih.gov/nuccore/CP001157 on 2020-05-06.

### Command Line

Entrez Direct (EDirect) has command line tools to interact with NCBI. To install these tools, run one of these lines in your command line. 
If having trouble, please reference [this guide](https://www.ncbi.nlm.nih.gov/books/NBK179288/). 
```bash
sh -c "$(curl -fsSL ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
sh -c "$(wget -q ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"
```

We can then retrieve the FASTA files using the following commands: 
```bash
efetch -db nuccore -id CP001189 -format fasta > CP001189.1.fasta
efetch -db nuccore -id CP001189 -format fasta -revcomp > CP001189.1_revcomp.fasta
efetch -db nuccore -id CP001190 -format fasta > CP001190.1.fasta
efetch -db nuccore -id CP001190 -format fasta -revcomp > CP001190.1_revcomp.fasta
```
> At the time of writing, this was not possible to do for GFF3 files.

## Part 1

The next part is a series of commands that search the FASTA sequence for TA sites and records their location. 
If doing this on your personal computer, open your command line app (on windows you can search "Command Prompt"). 
Activate bash by typing `bash` and hitting enter. 
If you get an error, google using bash on windows and follow the most current instructions.

To paste, use `CTRL + SHIFT + V`.

Make a gene_list.txt and run the shell script to find the TA sites for each "chromosome". 
```bash
for i in *.gff3;
do
  FILE=$(echo ${i} | sed "s/\.gff3//")
  python3 parse_gff3_file.py ${i} > ${FILE}_gene-list.txt
  ./analyze_TA-sites.sh ${FILE}
done
``` 

Add a column to each 
```bash
cat *nonpermissive.txt > TA-sites_nonpermissive.txt
cat *TA_flanks.txt > TA_flanks.txt
python3 find_homologous_TA_sites.py TA_flanks.txt > TA-sites_homologous.txt
```

Count the lines in each file to get the number of each type of TA site. 
Save this output to site_cnts.txt. 
```bash
wc -l *TA-sites*.txt | sed \$d > site_cnts.txt
```

Count the total number of genes in the organism. 
Save to gene_cnts.txt
```bash
cat *gene-list.txt | cut -f2 | uniq -c | wc -l | sed 's/$/ Total Genes/' > gene_cnts.txt
```

Count the total number of removed genes (i.e. those that have no useable TA sites). 
Append this result to gene_cnts.txt with `>>` instead of `>`.
```bash
cat *removed_genes.txt | cut -f2 | uniq -c | wc -l | sed 's/$/ Removed Genes/' >> gene_cnts.txt
```

Combine the annotated TA sites and fasta files for use in **Part 2**. 
```bash
cat *unfiltered.txt > TA-sites_unfiltered.txt
cat CP001189.1.fasta CP001190.1.fasta > gdj.fasta
```

## Part 2

Part two converts the sequencing data into counts.

- Input: demultiplexed fastq files from a single sequencing run
- Output: count matrix of TA insertions with gene labels

To get started, connect to the university's internet either through eduroam or VPN, then open PuTTY and login with your username and password. 
Once connected, we can organize our project folder. 
```bash
mkdir -p gdj-tnseq/ {raw-data_fastq,genome-info,scripts,final-counts_tsv,processed-seqs}
```

Then we can move into the data release folder and copy the sequencing files. 
After copying, we can move back to our home directory with `cd ~`.
```bash
cd /home/barneyb/data_release/umgc/novaseq/210210_A00223_0498_BH333KDRXY/Barney_Project_024
cp -v *R1_001.fastq.gz ~/gdj-tnseq/raw-data_fastq
cd ~
```

Next, we need to upload two files from **Part 1** to `genome-info/` and three scripts into the `scripts/` folder.
I use WinSCP to do this. 
The files we need are: 

- `gdj.fasta`
- `TA-sites_unfiltered.txt`
- `count_TAs.py`
- `get-counts.sh`
- `run_get-counts.sh`

To make it easier to get started, I uploaded my `sw` (software) folder to improve the reproducability of this code. 
To use, upload to your home directory using WinSCP then decompress by using the following command in PuTTY.
```bash
tar -zxvf sw.tar.gz 
```

When finished, your project folder should like this:
```bash
username/
├── gdj-tnseq/
│   ├── final-counts_tsv/
│   ├── genome-info/
│   │   ├── gdj.fasta
│   │   └── TA-sites_unfiltered.txt
│   ├── processed-seqs/
│   ├── raw-data_fastq/
│   └── scripts/
│       ├── count_TAs.py
│       ├── get-counts.sh
│       └── run_get-counts.sh
└── sw/
    ├── bioawk/
    └── bowtie-1.2.3/
```

Now we are ready to connect to the HPC system.
```bash
ssh mangi
```

Once connected to mangi, we need to build our reference files for doing alignments later on. 
Use the bowtie-build command to build the reference genome. 
Bowtie wasn't available as a module, so I installed it in my local folder.
```bash
cd ~/gdj-tnseq/genome-info
~/sw/bowtie-1.2.3/bowtie-build gdj.fasta gdj
```

From here we submit a job or do an interactive session.
We will be submitting a job, but I included the command for an interactive session just in case.

In an interactive sesssion, we enter the commands directly at the command line. 
To start an interactive session, we can enter the following.
```bash
qsub -I -l nodes=1:ppn=4,walltime=8:00:00,mem=8gb -q interactive
```
Type `exit` to leave an interactive session.

To submit a job, we will need scripts with the right permissions.
```bash
cd ~/gdj-tnseq/scripts/
chmod u+r+x *.sh
chmod u+r+x *.py
```

Once a job is submitted, it will wait in a queue until resources become available.
To submit our job, we simply type the command below. 
Before doing this though, make sure all your file paths are correct. 
```bash
cd ~/gdj-tnseq/scripts/
sbatch  run_get-counts.sh
```

**Both Bioawk and Bowtie are required for this script.
I installed these in my profile in a `~/sw/` directory. 
Before beginning, you will need to get these set up in your home directory.**

Some helpful commands: 
- `squeue -u username` to check the status of the job you submitted
- `scancel jobIDnumber` to cancel a job
- `groupquota` to make sure our group has enough storage to proceed

> Note: this script is not optimized ([might do this later](https://pages.github.umn.edu/dunn0404/job-submission-and-scheduling/16-resources/index.html))

Another thing I did was combine the alignment output for each file into one master file. 
```bash
head -n -0 *alignment-info.txt > _alignment-info_ALL.txt
```

Now we can download the filtered counts with WinSCP to our personal computer. 
These counts will be used to calculate fitness in R. 

## Part 3

This part is largely up to you and how you wish to calculate fitness as there are many methods used in the literature. 
We chose to use the van Opijnen method as the final values held a practical definition and are more digestable for those unfamiliar to Tn-seq analysis. 

The scripts I wrote are provided in the `02_calc-fitness` folder. 
Each script is written to be opened in R Studio and executed one at a time.
I preferred to not write these files as executables so that I could view and work with the data at each step, but now that the process has been determined, we could convert them. 
Each script is named to indicate the order in which they should be executed. 

## Part 4

This part was something we did for the Avn Tn-seq manuscript to compare essential genes between Avn and P.
I also did this for gdj and avn to assist in gene determination since the gdj genome is not as well annotated as avn.
To establish homologs between the two organisms, we used the desktop version of Blast. 
The first step is to download Blast for your computer: https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html.
The results can then be imported into R for further processing. 

```bash
# install miniconda
bash Miniconda3-latest-Linux-x86_64.sh
export PATH=~/miniconda3/condabin:$PATH

# Create environment 

# install blast
conda config --add channels bioconda
conda install blast

# get tabular gene list
python3 parse_gff3_file.py Azotobacter_vinelandii_DJ.gff3 > avn_gene_list.txt

# extract only essential protein seqs from list of all protein seqs
seqkit grep -n -r -f avn_essential.txt CP001157_protein-seqs.txt > avn_essential-seqs.fasta

# make ref database
makeblastdb -in CP001157_protein-seqs.fasta -dbtype prot

# blast protein seqs against eachother
blastp -query gdj_protein-seqs.fasta -db CP001157_protein-seqs.fasta -out gdj-avn_matches.tsv -outfmt 6

# header for the output data
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
```
