library(tidyverse)
library(here)
library(pheatmap)

parse_ariba_file <- function(filepath) {
  read_csv(filepath, show_col_types = FALSE) %>%
    mutate(name = gsub("./(redo_failed_)?ariba_amr_results/|/report.tsv", "", name)) %>%
    select(contains("name") |
      contains("match") |
      contains("assembled") |
      contains("ref_seq")) %>%
    pivot_longer(-name,
      names_to = c("cluster_name", ".value"),
      names_sep = "\\."
    ) %>%
    filter(match == "yes") %>%
    select(name, ref_seq) %>%
    mutate(ref_seq = sapply(strsplit(ref_seq, split = "__"), function(x) {
      x[3]
    })) -> dfout
  return(dfout)
}

ariba_dir <- "data-raw/ecoli-genomics-paper/horesh_full_diversity/ariba_amr/"

ariba_files <- paste0(ariba_dir,list.files(ariba_dir))
           
map_df(ariba_files, parse_ariba_file) -> btESBL_ecoli_horesh_amr

usethis::use_data(btESBL_ecoli_horesh_amr, overwrite = TRUE)
