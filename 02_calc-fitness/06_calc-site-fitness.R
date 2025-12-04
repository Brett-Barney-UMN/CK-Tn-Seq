
# Calculates the fitness for each site (minus non-unique or non-permissive)

library(tidyverse)


# Load data ---------------------------

read_tsv(here::here("./data/processed/t1_site-avgs.tsv")) -> t1_site_avgs

read_tsv(here::here("./data/processed/t2_site-cnts.tsv")) -> t2_site_cnts

read_tsv(here::here("./data/processed/expansion-factors.tsv")) -> exp_factors


# Create fitness equation ---------------------------

fitness <- function(t1, t2, d) {
  result <- log(t2 * (d / t1)) / log((1 - t2) * (d / (1 - t1)))
  return(result)
}


# Calculate site fitness ---------------------------

# change column name for t2 counts
names(t2_site_cnts)[names(t2_site_cnts) == "site_cnt"] <- "t2_site_cnt"

# combine data 
t2_site_cnts %>%
  left_join(t1_site_avgs) %>% 
  left_join(exp_factors %>% select(-od600, -org), by = c("media", "rep")) %>%
  group_by(media, rep) %>% 
  mutate(
    t1_site_p = t1_site_cnt / sum(t1_site_cnt), 
    t2_site_p = t2_site_cnt / sum(t2_site_cnt), 
    w = fitness(t1 = t1_site_p, t2 = t2_site_p, d = d)
    ) %>% 
  ungroup() %>% 
  select(-d, -starts_with("t2"), -t1_site_p, -description) -> 
  site_fit 

# pivot to wide format for nitrogen manuscript
site_fit %>% 
  mutate(w = round(w, digits = 3)) %>% # round numbers to three decimal places
  pivot_wider(
    names_from = c(media, rep),  
    values_from = w
  ) -> 
  site_fit_wide


# Write to file ---------------------------

site_fit %>% write_tsv(here::here("./data/processed/fitness_per-site.tsv"))

site_fit_wide %>% write_tsv(here::here("./data/processed/fitness_per-site_wide.tsv"))
