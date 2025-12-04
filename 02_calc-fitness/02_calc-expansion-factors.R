
# Calculates the expansion factor for each sample from OD600 data 

library(tidyverse)

# get OD data 
ods <- read_tsv(file = "./data/raw/od600.txt") 

# clean up column names
colnames(ods) <- colnames(ods) %>% str_replace_all("\\s", "_") %>% tolower()

# calculate expansion factor for each sample

od_t1 <- 0.025 # starting OD for each culture

ods %>% 
  select(id, od600) %>%
  separate(id, into = c("org", "media", "rep"), sep = "_") %>% 
  filter(media != "Initial") %>% 
  mutate(d = od600/od_t1) -> 
  exp_factors

exp_factors %>% write_tsv(here::here("./data/processed/expansion-factors.tsv"))
