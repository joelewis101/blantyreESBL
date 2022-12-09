library(tidyverse)
library(here)
library(pheatmap)

parse_ariba_file <- function(filepath, show_col_types = FALSE) {
  read_csv(filepath, show_col_types = show_col_types) %>%
    mutate(
      name =
        gsub(
          "./(redo_failed_)?(ariba_amr_results)?/|^\\./|/report.tsv",
          "",
          name
        )
    ) %>%
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

ariba_files <- here(paste0(ariba_dir, list.files(ariba_dir)))

map_df(ariba_files, parse_ariba_file) -> btESBL_ecoli_horesh_amr


btESBL_ecoli_horesh_amr %>%
  mutate(name = gsub("_contigs_.*$", "", name)) %>%
  left_join(
    read_tsv(
      here(
        "data-raw/ecoli-genomics-paper/horesh_full_diversity/FINAL_METADATA_CLEANED.tsv"
      )
    ) %>%
      separate(Reads_Location,
        into = paste0(
          1:4,
          "reads"
        ),
        sep = ","
      ) %>%
      select(ID, contains("reads")) %>%
      pivot_longer(-ID) %>%
      group_by(ID) %>%
      arrange(ID, value) %>%
      slice(1) %>%
      mutate(value = gsub("^/.*/|(_1)?\\.fastq(\\.gz)?", "", value)) %>%
      select(-name),
    by = c("name" = "value")
  ) %>%
  rename(read_id = name) -> btESBL_ecoli_horesh_amr


usethis::use_data(btESBL_ecoli_horesh_amr, overwrite = TRUE)
