
# removes non-unique sites and checks non-permissive sites

library(tidyverse)


# Load data ---------------------------

# counts
site_cnts <- read_tsv("./data/processed/site-counts.tsv")

# define header for nonpermissive and nonunique sites
fps_header <- c("seq", "t_pos", "chrm")

read_tsv(
  file = "./data/processed/TA-sites_homologous.txt", 
  col_names = fps_header) -> 
  nonuniq

read_tsv(
  file = "./data/processed/TA-sites_nonpermissive.txt", 
  col_names = fps_header) %>% 
  mutate(type = "nonperm") -> 
  nonperm


# Check insertion rates for non-permissive sites ---------------------------

# should be 58303 total TA sites
# here are numbers from an old email I sent to Erin and Brett
57574 + 729

# number of non-unique sites that would be removed
# some sites are duplicates because the fwd and rev flank are both non-unique
nonuniq %>% 
  select(-seq) %>% 
  unique() %>% 
  nrow()

site_cnts %>% 
  # remove non-unique sites
  anti_join(nonuniq, by = c("chrm", "t_pos")) %>%
  # label non-permissive and permissive sites
  left_join(nonperm %>% select(-seq), by = c("chrm", "t_pos")) %>%  
  filter(media == "Initial") -> 
  site_cnts_uniq

# how many sites are in each group
site_cnts_uniq %>% count(rep, type)

# how many non-zero sites are in each group
site_cnts_uniq %>% filter(site_cnt > 0) %>% count(media, rep, type)

# calculate the mean with standard error for the site cnt 
# and also for the non-zero site cnts
site_cnts_uniq %>% 
  group_by(type) %>% 
  summarise(
    avg_site_cnt = mean_se(site_cnt), 
    no_zeros = mean_se(site_cnt[site_cnt > 0])
  ) %>% 
  rename() -> 
  mean_se_results 

# get the summary stats for each permissive and non-permissive sites
site_cnts_uniq %>% 
  select(site_cnt, type) %>% 
  split(.$type) %>% 
  map(summary)


# Plot results ---------------------------

# set plot theme
theme_set(theme_bw(base_size = 9))
theme_update(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), 
  panel.spacing.y = unit(0, "lines"), # no vertical space between facets
  strip.text.y = element_text(angle = 0), # facet labels
  legend.position = "none",
  plot.margin = unit(c(2, 2, 2, 2), "mm"), 
  strip.background = element_blank())


# plot mean and standard error
site_cnts_uniq %>% 
  ggplot(aes(x = grp, y = site_cnt)) + 
  geom_jitter(col = "lightblue", alpha = 0.6, width = 0.1) +
  theme_bw() + 
  facet_grid(. ~ rep) + 
  labs(x = NULL) +
  theme(axis.text.x = element_blank()) -> 
  p1

p1

# plot mean and standard error for non-zero sites
site_cnts_uniq %>% 
  filter(site_cnt > 0) %>% 
  ggplot(aes(x = grp, y = site_cnt)) + 
  geom_jitter(col = "lightblue", alpha = 0.6, width = 0.1) +
  theme_bw() + 
  facet_grid(. ~ rep) -> 
  p2

p2

# combine plots 
(p1 + p2) + plot_layout(ncol = 1) + plot_annotation(tag_levels = 'A') -> g

ggsave("./docs/figs/nonpermissive-sites.tiff", 
       plot = g, dpi = 400, width = 168, height = 125, units = "mm")


# Write to file ---------------------------

site_cnts %>% 
  anti_join(nonuniq, by = c("chrm", "t_pos")) %>%
  anti_join(nonperm, by = c("chrm", "t_pos")) %>%
  write_tsv("./data/processed/site-counts_filtered.tsv")
