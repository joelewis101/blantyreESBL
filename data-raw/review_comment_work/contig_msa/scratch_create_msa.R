# plot multiple sequence alignment of genes

# first - decide on what clusters to plot
library(here)
library(tidyverse)
library(blantyreESBL)

btESBL_contigclusters %>%
  group_by(clstr_name) %>%
  tally() %>%
  arrange(-n)

btESBL_contigclusters %>%
  filter(grepl("CTX_M_27", clstr_name)) %>%
  ggplot(aes(length)) +
  geom_histogram(bins = 100) +
  facet_wrap(~clstr_name, scales = "free")


btESBL_contigclusters %>%
  filter(clstr_name == "CTX_M_27.1") %>%
  group_by(length) %>%
  tally()

btESBL_contigclusters %>%
  group_by(clstr_name) %>%
  mutate(n = n())  %>%
  filter(grepl("CTX_M_15", clstr_name), n > 3) %>%
  mutate(clstr_name = as.factor(clstr_name),
         clstr_nme = fct_reorder(clstr_name, n, .desc = TRUE)) %>%
  ggplot(aes(length)) +
  geom_histogram(bins = 100) +
  facet_wrap(clstr_nme ~ ., scales = "free")


# ok plot top 3 largest
btESBL_contigclusters %>%
  group_by(clstr_name) %>%
  tally() %>%
  arrange(-n) %>%
  slice(1:10) %>%
  pull(clstr_name) -> clusters_for_msa

slice_and_save_clusters <- function(df, cluster, path = "") {
  df %>%
    filter(clstr_name == {{ cluster }},
           clstr_rep == 0) %>%
    dplyr::pull(id) -> contigs_to_save
  write_lines(
    contigs_to_save,
    paste0(path, cluster, ".txt")
  )
  df %>%
    filter(clstr_name == {{ cluster }},
           clstr_rep == 1) %>%
    dplyr::pull(id) -> contigs_rep_to_save

  write_lines(
    contigs_rep_to_save,
    paste0(path, cluster, "_rep.txt")

  )
}



lapply(clusters_for_msa, slice_and_save_clusters, df = btESBL_contigclusters,
       "data-raw/review_comment_work/contig_msa/")

# what about 0.99 seq identity

source(
  here(paste0(
      "data-raw/review_comment_work/",
      "contig_sens_ax/scratch_load_contig_sens_ax_data.R" )
      )
  )

df %>%
  filter(len_diff_cutoff == 0.8, ident_cutoff == 0.95) %>%
  mutate(clstr_name = paste0(gene, ".", cluster)) %>%
  group_by(clstr_name) %>%
  tally() %>%
  arrange(-n)

df %>%
  filter(len_diff_cutoff == 0.8, ident_cutoff == 0.95) %>%
  mutate(clstr_name = paste0(gene, ".", cluster)) %>%
  group_by(clstr_name) %>%
  tally() %>%
  arrange(-n) %>%
  slice(1:3) %>%
  pull(clstr_name) -> clusters_for_msa


lapply(clusters_for_msa,
  slice_and_save_clusters,
  df = df %>%
    filter(len_diff_cutoff == 0.8, ident_cutoff == 0.95) %>%
    mutate(
      clstr_name = paste0(gene, ".", cluster),
      id = paste0(".", contig)),
  path = "data-raw/review_comment_work/contig_msa/len_cutoff0.8/"
)


#### for minimap

do.call(
  bind_rows,
  map(
    list.files(
      here("data-raw/contig_clusters/"),
      pattern = "tsv",
      full.names = TRUE
    ),
    ~ read_tsv(.x,
               col_types = "ccccccc",
               id = "file"
    )
  )
) -> df_orig

df_orig %>%
  mutate(
    gene = str_extract(
    file,
    "(?<=contigs-).+(?=\\.clust)"
  ),
  clstr_name = paste0(gene, ".", clstr)) -> df_orig



df_orig %>%
  group_by(clstr_name) %>%
  mutate(n = n()) %>%
  arrange(desc(n), clstr_name)

df_orig %>%
  filter(clstr_name == "CTX_M_27.1", clstr_rep == 1 ) %>%
  pull(id) -> ctxm27.1_cluster_rep.txt

df_orig %>%
  filter(clstr_name == "CTX_M_27.1", clstr_rep == 0 ) %>%
  pull(id) -> ctxm27.1_others.txt

write_lines(
  ctxm27.1_cluster_rep.txt,
  "data-raw/review_comment_work/mapping/ctxm27.1_cluster_rep.txt"
)

write_lines(
  ctxm27.1_others.txt,
  "data-raw/review_comment_work/mapping/ctxm27.1_others.txt"
)

# ARE these the same things?

df_orig %>%
  filter(clstr_name == "CTX_M_27.1") %>%
  pull(id) -> ctxm27.1.txt

btESBL_contigclusters %>%
  filter(clstr_name == "CTX_M_27.1") %>%
  pull(id) -> b

