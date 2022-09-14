library(tidyverse)
library(here)

st410_snpdists <-
  read_tsv(
    here("data-raw/ecoli-genomics-paper/st410-new/st410-snpdists.tsv")
  )

names(st410_snpdists)[1] <- "name1"

st410_malawi_lanes <- read_lines(
  here("data-raw/ecoli-genomics-paper/st410-new/st410-malawi-clade-samples.txt")
)

st410_snpdists %>%
  pivot_longer(-name1) %>%
  filter(name1 %in% st410_malawi_lanes,
         name %in% st410_malawi_lanes,
         name1 != name) %>%
  pull(value) |>
quantile(c(0,0.25,0.5,0.75,1))



st167_snpdists <-
  read_tsv(
    here("data-raw/ecoli-genomics-paper/st167-new/st167-snpdists.tsv")
  )

names(st167_snpdists)[1] <- "name1"

st167_malawi_lanes <- read_lines(
  here("data-raw/ecoli-genomics-paper/st167-new/st167-malawi-clade-samples.txt")
)

st167_snpdists %>%
  pivot_longer(-name1) %>%
  filter(name1 %in% st167_malawi_lanes,
         name %in% st167_malawi_lanes,
         name1 != name) %>%
  pull(value) |>
quantile(c(0,0.25,0.5,0.75,1))
