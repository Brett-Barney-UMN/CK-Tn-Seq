
# Combines all the filtered count files from the MSI pipeline into one file 
# and combines the fwd and rev counts together. Filtered counts have already had the 
# non-permissive sites and non-unique sites removed. 

library(tidyverse)


# Load data ---------------------------

list.files(here::here("./data/processed/unfiltered-counts/"), pattern = "txt", full.names = TRUE) -> 
  fps

# define header for filtered counts
fps_header <- c("chrm", "t_pos", "locus_tag", "fwd_count", "rev_count")

# read in the filtered counts as a nested tibble
tibble(
  media = basename(fps) %>% strsplit(., "_") %>% lapply(., `[[`, 2) %>% unlist(),
  rep = basename(fps) %>% strsplit(., "_") %>% lapply(., `[[`, 3) %>% unlist(),
  cnts = map(.x = fps, .f = read_tsv, col_names = fps_header)) -> 
  nested_data

# Process data ---------------------------

# combine fwd and rev counts for each TA site
nested_data %>% 
  unnest(cols = c(cnts)) %>%
  mutate(site_cnt = fwd_count + rev_count) %>%
  select(c(-fwd_count, -rev_count)) ->
  unnested_data


# Write to file ---------------------------

unnested_data %>% write_tsv(here::here("./data/processed/site-counts.tsv")) 
