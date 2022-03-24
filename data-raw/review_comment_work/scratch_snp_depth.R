library(tidyverse)
library(here)
library(blantyreESBL)

bind_rows(
read_tsv("~/Documents/PhD/Manuscripts/esbl_carriage/manuscript/nature_microbiology/snp_depths/snp_depths_all_esco.tsv") %>%
  mutate(species = "E. coli"),
read_tsv("~/Documents/PhD/Manuscripts/esbl_carriage/manuscript/nature_microbiology/snp_depths/snp_depths_all_kleb.tsv") %>%
  mutate(species = "K. pneumoniae complex")
) -> snp_depths

snp_depths  %>%
  group_by(species) %>%
  summarise(mean_dep = mean(snp_depth, na.rm = TRUE),
            sd = sd(snp_depth, na.rm = TRUE))
#
# snp_depths %>%
#   group_by(position) %>%
# mutate(mean_dep = mean(snp_depth, na.rm = TRUE),
#             sd = sd(snp_depth, na.rm = TRUE)) %>%
#   ggplot(aes(position, mean_dep)) + geom_point()

snp_depths %>%
ggplot(aes(snp_depth, fill = species)) +
  geom_histogram(binwidth = 1, position = "dodge" )  +
  coord_cartesian(xlim=c(0,100)) + theme_bw()

bind_rows(
  read_tsv(here("data-raw/review_comment_work/snp_depth/core_esco.txt")) %>%
  mutate(species = "E. coli"),
  read_tsv(here("data-raw/review_comment_work/snp_depth/core_kleb.txt")) %>%
  mutate(species = "K. pneumoniae complex")) ->
  core_all

core_all %>%
  ungroup() %>%
  summarise(med = median(LOWCOV),
         lq = quantile(LOWCOV, 0.25))

core_all %>%
  mutate(prop_ambiguous = LOWCOV/ALIGNED) %>%
  ggplot(aes(prop_ambiguous, fill = species)) +
  geom_histogram(bins = 1000, position = "dodge") +
  coord_cartesian(xlim=c(0,0.1)) + theme_bw()

core_all %>%
  mutate(prop_ambiguous = LOWCOV/ALIGNED) %>%
  group_by(species) %>%
  summarise(median = median(prop_ambiguous),
            lci = quantile(prop_ambiguous, 0.25),
            uci = quantile(prop_ambiguous, 0.75))

  summarise(mean = mean(prop_ambiguous),
            sd = sd(prop_ambiguous))


core_all %>%
  mutate(prop_ambiguous = LOWCOV/VARIANT) %>%
  ggplot(aes(prop_ambiguous, fill = species)) +
  geom_histogram(bins = 1000, position = "dodge") +
  coord_cartesian(xlim=c(0,100)) + theme_bw()

snp_dist_kleb_a <- read_tsv(
  here("data-raw/review_comment_work/snp_depth/snp_dists_kleb_alli_opt_a.tsv")
) %>%
  rename(sample_a = `snp-dists 0.6.2`)

snp_dist_kleb_a %>%
  pivot_longer(-sample_a,
               names_to = "sample_b") -> snp_dist_kleb_a

bind_rows(
  snp_dist_kleb_a %>%
    mutate(type = "sens ax"),
  btESBL_snpdists_kleb %>%
    rename(sample_a = sample) %>%
    pivot_longer(-sample_a,
                 names_to = "sample_b") %>%
  mutate(type = "original")
  ) %>%
  ggplot(aes(value, group = type, fill = type)) +
    geom_histogram(bins = 100,position = "dodge") +
  coord_cartesian(xlim = c(0,100000))

