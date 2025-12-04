
# Calculates the fitness for each gene (only sites between 05-90% used)

library(tidyverse)


# Load data ---------------------------

read_tsv(here::here("./data/processed/t1_gene-avgs.tsv")) -> t1_gene_avgs

read_tsv(here::here("./data/processed/t2_gene-cnts.tsv")) -> t2_gene_cnts

read_tsv(here::here("./data/processed/expansion-factors.tsv")) -> exp_factors


# Create fitness equation ---------------------------

fitness <- function(t1, t2, d) {
  result <- log(t2 * (d / t1)) / log((1 - t2) * (d / (1 - t1)))
  return(result)
}


# Calculate gene fitness ---------------------------

# combine data 
t2_gene_cnts %>% 
  left_join(t1_gene_avgs) %>% 
  left_join(exp_factors %>% select(-org, -od600), by = c("media", "rep")) %>%
  group_by(media, rep) %>% 
  mutate(
    t1_gene_p = t1_gene_cnt / sum(t1_gene_cnt), 
    t2_gene_p = t2_gene_cnt / sum(t2_gene_cnt), 
    w = fitness(t1 = t1_gene_p, t2 = t2_gene_p, d = d)
    ) %>% 
  ungroup() %>% 
  select(-d, -starts_with("t2"), -t1_gene_p) -> 
  gene_fit


# Assign categories ---------------------------

gene_fit %>% 
  group_by(media, locus_tag) %>%
  summarise(
    w_avg = mean(w), 
    .groups = "drop"
    ) %>%
  mutate(cat = case_when( 
    # anything not fitting in one of the categories gets NA
    w_avg > 1.10 & w_avg != Inf   ~ "growth_promoting",
    w_avg <= 1.10 & w_avg > 0.90  ~ "none",
    w_avg <= 0.90 & w_avg > 0.55  ~ "mod_defect",
    w_avg <= 0.55 & w_avg >= 0    ~ "lrg_defect",
    w_avg < 0                     ~ "lrg_defect", 
    TRUE                          ~ "no_data")
    ) -> 
  gene_cats


# Write to file ---------------------------

gene_fit %>%
  # add fitness categories
  left_join(gene_cats %>% select(locus_tag, media, cat), by = c("locus_tag", "media")) %>%
  arrange(locus_tag) %>% 
  # reorder columns
  select(locus_tag, name, description, t1_gene_cnt, t1_n_ta, t1_n_hit, everything()) -> 
  fit_output

fit_output %>% write_tsv(here::here("./data/processed/fitness_per-gene.tsv"))


# Wide version for providing others

fit_output %>% 
  select(-cat) %>% 
  mutate(
    w = round(w, 3), 
    w = as.character(w)
  ) %>% 
  pivot_wider(
    names_from = c(media, rep),  
    values_from = w, 
  ) -> 
  fit_output_wide

fit_output_wide %>% write_tsv(here::here("./data/processed/fitness_per-gene_wide.tsv"))

