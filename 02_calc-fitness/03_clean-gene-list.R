
# Clean the gene list for avn by correcting punctuation and updating gene names
# per the literature.


library(tidyverse)


# Load data ---------------------------

# create header for gene list
file_header <- c("chrm", "locus_tag", "name", "description", "start_pos", "end_pos", "strand")

# load gene list
read_tsv(
  file = "./data/processed/gene-list_gdj.txt", 
  col_names = file_header
  ) %>% 
  # correct punctuation
  mutate(
    description = gsub("%2C", ",", description),
    description = gsub("%3B", ";", description)
    ) -> 
  gene_list

# Process data ---------------------------

# how many duplicate entries do we have?
gene_list$locus_tag[duplicated(gene_list$locus_tag)] %>% length()

# combine duplicate entries so that there is only 
# one entry in the gene list per locus_tag
gene_list %>% 
  group_by(locus_tag) %>% 
  mutate(
    start_pos = min(start_pos), 
    end_pos = max(end_pos)
    ) %>% 
  ungroup() %>% 
  distinct() -> 
  gene_list_cln

# check to see if the duplicate entries are removed
gene_list_cln$locus_tag[duplicated(gene_list_cln$locus_tag)] %>% length()


# Save output ---------------------------

# write to file
gene_list_cln %>% write_tsv(here::here("./data/processed/gene-list_gdj_cln.tsv"))
