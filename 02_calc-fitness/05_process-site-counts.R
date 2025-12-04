
# calculates the location of each TA site in the gene and outputs the various files 
# needed for downstream analysis 

library(tidyverse)


# Load data ---------------------------

site_cnts <- read_tsv(file = "./data/processed/site-counts_filtered.tsv") 

gene_list <- read_tsv(here::here("./data/processed/gene-list_gdj_cln.tsv"))


# Determine site locations ---------------------------

site_cnts %>%
  left_join(gene_list %>% select(chrm, locus_tag, name, description, start_pos, end_pos, strand)) %>% 
  mutate(
    ins_location = round((t_pos - start_pos) / (end_pos - start_pos), digits = 3),
    ins_location = case_when( # correct insertion locations for genes on the rev strand
      strand == "-" ~ 1 - ins_location, 
      strand == "+" ~ ins_location
    )) %>%
  select(-start_pos, -end_pos, -strand) %>% 
  select(chrm, locus_tag, name, description, t_pos, ins_location, media, rep, site_cnt) -> 
  site_locs


# Average t1 site counts ---------------------------

site_locs %>% 
  filter(media == "Initial") %>% 
  group_by(rep) %>% 
  summarise(total_reads = sum(site_cnt)) ->
  total_reads

site_locs %>% 
  filter(media == "Initial") %>% 
  group_by(rep) %>% 
  mutate(cpm = site_cnt / sum(site_cnt) * mean(total_reads$total_reads)) %>%  # counts per million
  ungroup() %>% 
  group_by(chrm, locus_tag, name, description, t_pos, ins_location) %>% 
  summarise(
    t1_site_cnt = mean(cpm)
  ) %>% 
  ungroup() ->
  t1_site_avgs
  

# Create output files ---------------------------

site_locs %>% 
  filter(media != "Initial") %>% 
  select(-description) -> 
  t2_site_cnts

t1_site_avgs %>% 
  filter(ins_location < 0.9, ins_location > 0.05) %>%  
  group_by(locus_tag, name, description) %>% 
  summarise(
    t1_n_ta     = n(),
    t1_n_hit    = sum(t1_site_cnt > 0),
    t1_gene_cnt = sum(t1_site_cnt),
    .groups = "drop"
  ) -> 
  t1_gene_avgs

site_locs %>% 
  filter(media != "Initial") %>% 
  filter(ins_location < 0.9, ins_location > 0.05) %>% 
  group_by(media, rep, locus_tag) %>% 
  summarise(
    t2_gene_cnt = sum(site_cnt),
    .groups = "drop"
  ) -> 
  t2_gene_cnts


# Write to file ---------------------------

# t1 site averages and t2 sites output for site fitness calculations

t1_site_avgs %>% write_tsv(here::here("./data/processed/t1_site-avgs.tsv"))
t2_site_cnts %>% write_tsv(here::here("./data/processed/t2_site-cnts.tsv"))


# t1 and t2 gene output for essential gene analysis (t1) and gene fitness calculations (t1 and t2)
# data has tip and tail end sites removed

t1_gene_avgs %>% write_tsv(here::here("./data/processed/t1_gene-avgs.tsv"))
t2_gene_cnts %>% write_tsv(here::here("./data/processed/t2_gene-cnts.tsv"))


# previous code to look at different filtering
# gene_p_t2_0090 = sum(site_cnt[ins_location < 0.9 & ins_location > 0.00]),
# gene_p_t2_0590 = sum(site_cnt[ins_location < 0.9 & ins_location > 0.05]),
# gene_p_t2_1090 = sum(site_cnt[ins_location < 0.9 & ins_location > 0.10]),
