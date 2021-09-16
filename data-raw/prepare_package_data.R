# Load raw data and save it as .rda files for blantyreESBL package

# the load_phd_data scripts are from the Thesis repo https://github.com/joelewis101/thesis
#library(tidyverse)
library("phytools")

source("/Users/joelewis/Documents/PhD/Thesis/bookdown/final_cleaning_scripts/load_PhD_data.R")

source("/Users/joelewis/Documents/PhD/Thesis/bookdown/final_cleaning_scripts/load_and_clean_lims.R")

library("here")
library("ape")

source(here("data-raw/load_metadata_fn.R"))

# baseline characteristics -------------------------------------------

enroll %>%
  filter(arm != 4) %>%
  mutate(
    art_time =  as.numeric(
      (data_date - hivartstart) / (365.25 /  12)
      ),
    hivart = if_else(hivart != "5A", "ART regimen: other", hivart),
    animalskept = if_else(animalskept == "Other", "Elsewhere", animalskept),
    arm = as.character(arm),
    enroll_date = data_date,
  ) %>% select(pid,
    arm,
    enroll_date,
    calc_age,
    ptsex,
    hivstatus,
    recieved_prehosp_ab,
    pmhxrechospital,
    tbstatus,
    tbongoing,
    hivonart,
    hivart,
    art_time,
    hivcpt,
    tobaccoyn,
    alcoholyn,
    highestedu,
    job,
    housholdadultsno,
    householdchildno,
    toilet ,
    watersource,
    watertreated,
    electricityyn,
    fuel,
    keepanim,
    keep.poultry,
    keep.goats,
    keep.dogs,
    keep.cattle,
    keep.sheep,
    keep.mules
  ) -> btESBL_participants


use_data(btESBL_participants, overwrite = TRUE)

### longitudinal exposure (covariate) data ---------------------------------

read_csv(
  "~/Documents/PhD/Thesis/bookdown/data/longit_covariate_data.csv") ->
  btESBL_exposures

# add in cotrim

btESBL_exposures %>%
  left_join(select(enroll, pid, hivcpt) %>%
              unique()) %>%
  mutate(cotri = if_else(hivcpt == "Yes" & !is.na(hivcpt)
                           , 1, cotri)) %>%
  select(-hivcpt) -> btESBL_exposures

use_data(btESBL_exposures, overwrite = TRUE)

# ESBL carriage data -------------------------------------------------


lims_dates %>%
  rename("t" = "assess_type") ->
  btESBL_stoolESBL


use_data(btESBL_stoolESBL, overwrite = TRUE)

lims_orgs %>%
  select(lab_id,organism, ESBL) %>%
  # recode - negatives were  coded as all sorts
  mutate(ESBL = if_else(ESBL == "Positive",
                        "Positive",
                        "Negative")) %>%
  unique() %>%
  filter(!is.na(organism)) -> btESBL_stoolorgs

use_data(btESBL_stoolorgs, overwrite = TRUE)

# stan_data

read_csv("/Users/joelewis/Documents/PhD/Thesis/bookdown/data/stan_df.csv") ->
  btESBL_modeldata

btESBL_modeldata %>%
  select(-X1) -> btESBL_modeldata

use_data(btESBL_modeldata)

# stan models ----------------------------------------------

readRDS("/Users/joelewis/Documents/PhD/Thesis/bookdown/chapter_9/stan_models/model_2/stan_data_m2.rds") ->
  btESBL_model2data

use_data(btESBL_model2data, overwrite = TRUE)

readRDS("/Users/joelewis/Documents/PhD/Thesis/bookdown/chapter_9/stan_models/model_1/stan_data_m1.rds") ->
  btESBL_model1data

use_data(btESBL_model1data, overwrite = TRUE)

readRDS("/Users/joelewis/Documents/PhD/Thesis/bookdown/chapter_9/stan_models/model_2/stanfit_m2.rds") ->
  btESBL_model2posterior

use_data(btESBL_model2posterior, overwrite = TRUE)

readRDS("/Users/joelewis/Documents/PhD/Thesis/bookdown/chapter_9/stan_models/model_1/stanfit_m1.rds") ->
  btESBL_model1posterior

use_data(btESBL_model1posterior, overwrite = TRUE)

# simulations from posterior

btESBL_model2simulations <- read_csv("~/Documents/PhD/Thesis/bookdown/chapter_9/simulations2.csv")

btESBL_model2simulations %>%
  select(-c(pid, start_state, abx_cpt, tb_start))

use_data(btESBL_model2simulations, overwrite = TRUE)

# plasmid replicons

esco.plasm <- read_csv(here("data-raw/plasmids/esco-plasmidfinder-summ.csv"))
klebs.plasm <- read_csv("data-raw/plasmids/klebs-plasmidfinder-summary.csv")

dassim.klebs <- read_lines(here("data-raw/included_lanes/kleb_lanes_retained_following_qc.txt"))

bind_rows(

esco.plasm %>%
  pivot_longer(-name,
               names_to = c("cluster_name", ".value"),
               names_sep = "\\.") %>%
  mutate(name = gsub("\\./","", name),
         name = gsub("/report\\.tsv","",name)) %>%
  filter(match == "yes") %>%
  select(name, ref_seq) %>%
  mutate(species = "E. coli"),

klebs.plasm %>%
  mutate(name = gsub("\\./","", name),
         name = gsub("/report\\.tsv","",name)) %>%
  filter(name %in% dassim.klebs) %>%
  pivot_longer(-name,
               names_to = c("cluster_name", ".value"),
               names_sep = "\\.") %>%
  filter(match == "yes") %>%
  select(name, ref_seq) %>%
  mutate(species = "K. pneumoniae")
) %>%
  rename("sample" = "name") -> btESBL_plasmidreplicons


use_data(btESBL_plasmidreplicons, overwrite = TRUE)

# amr genes ---------------------------------------------------------


esco1 <- read_csv(here("data-raw/amr_genes/esco-amr-ariba-report-part1.csv"))
esco2 <- read_csv(here("data-raw/amr_genes/esco-amr-ariba-report-part2.csv"))
klebs <- read_csv(here("data-raw/amr_genes/kleb-amr-ariba-report.csv"))

bind_rows(
esco1 %>%
  select(contains("name") |
           contains("match") |
           contains("assembled") |
           contains("ref_seq")) %>%
  pivot_longer(-name,
               names_to = c("cluster_name", ".value"),
               names_sep = "\\.") %>%
  filter(match == "yes") %>%
  select(name, ref_seq) %>%
  mutate(ref_seq = sapply(strsplit(ref_seq, split = "__"), function(x)
    x[3]),
    species = "E. coli"),

esco2 %>%
  select(contains("name") |
           contains("match") |
           contains("assembled") |
           contains("ref_seq")) %>%
  pivot_longer(-name,
               names_to = c("cluster_name", ".value"),
               names_sep = "\\.") %>%
  filter(match == "yes") %>%
  select(name, ref_seq) %>%
  mutate(ref_seq = sapply(strsplit(ref_seq, split = "__"), function(x)
    x[3]),
    species = "E. coli"),

klebs %>%
  mutate(name = gsub("\\./","", name),
         name = gsub("/report\\.tsv","",name)) %>%
  filter(name %in% dassim.klebs) %>%
  select(contains("name") |
           contains("match") |
           contains("assembled") |
           contains("ref_seq")) %>%
  pivot_longer(-name,
               names_to = c("cluster_name", ".value"),
               names_sep = "\\.") %>%
  filter(match == "yes") %>%
  select(name, ref_seq) %>%
  mutate(ref_seq = sapply(strsplit(ref_seq, split = "__"), function(x)
    x[3]),
    species = "K. pneumoniae"
    )
) %>%
  mutate(name = gsub("\\./","", name),
         name = gsub("/report\\.tsv","",name),
         name = gsub("_1\\.fastq\\.gz", "", name)
  )  %>%
  rename("sample" = "name") -> btESBL_amrgenes

use_data(btESBL_amrgenes, overwrite = TRUE)

# core gene trees --------------------------------------------------------

read.tree(here("data-raw/core_gene_trees/esco_core_gene_tree.treefile")) ->
  btESBL_coregene_tree_esco

read.tree(here("data-raw/core_gene_trees/kleb_core_gene_tree.treefile")) ->
  btESBL_coregene_tree_kleb

btESBL_coregene_tree_esco <- midpoint.root(btESBL_coregene_tree_esco)
btESBL_coregene_tree_kleb <- midpoint.root(btESBL_coregene_tree_kleb)

use_data(btESBL_coregene_tree_esco, overwrite = TRUE)
use_data(btESBL_coregene_tree_kleb, overwrite = TRUE)

# sequence sample metadata ----------------------------------

# If there is no lane accession this will use sample accession
# Need to update once klebs go into the wild!

samp_metadata <- load_DASSIM3_metadata(
  location_of_phd_loading_script =
    "/Users/joelewis/Documents/PhD/Thesis/bookdown/final_cleaning_scripts/load_PhD_data.R",
  location_of_lims_loading_script =
    "/Users/joelewis/Documents/PhD/Thesis/bookdown/final_cleaning_scripts/load_and_clean_lims.R",
  lanes_to_include =
    "/Users/joelewis/Documents/Sanger/DASSIM3/data_raw/metadata/all_DASSIM_kleb_and_esco_lanes.txt",
  sanger_metadata_file <- "/Users/joelewis/Documents/Sanger/DASSIM3/data_raw/metadata/all_DASSIM_kleb_and_esco_metadata.csv",
  recent_dc_def = 120 )

outcome %>%
  group_by(pid) %>%
  arrange(pid,hospoutcome,hospoutcomedate) %>%
  slice(n = 1) -> outcome

samp_metadata %>%
  rename_with(~ tolower(gsub(" |\\.","_", .x))) %>%
  mutate(lane = gsub("#","_", lane)) %>%
  left_join(outcome,
            by = "pid") %>%
  select(lane, supplier_name,pid, arm, visit, data_date, enroll_date, assess_type,
         hosp_assoc, hospoutcomedate) -> btESBL_sequence_sample_metadata

left_join(
  btESBL_sequence_sample_metadata,
  bind_rows(
    read_csv(
      "~/Documents/PhD/Thesis/bookdown/chapter_7/global_tree/all_acc.csv"
    ),
    read_csv(
      "~/Documents/PhD/Manuscripts/20200208vanilla_genomics/carriage_genomics/data_raw/DASSIM3_accession.csv"
    )
  ) %>%
    mutate(
      `Lane name` = gsub("#", "_", `Lane name`),
      `Lane accession` = if_else(`Lane accession` == "not found",
                                 `Sample accession`,
                                 `Lane accession`)
    ) %>%
    transmute(lane = `Lane name`,
              accession = `Lane accession`),
  by = "lane"
)  %>%
  relocate(accession, before = everything()) ->
  btESBL_sequence_sample_metadata

# add poppunk clusters --------------------

bind_rows(
  read_csv(here(
    "data-raw/poppunk/strain_db_clustersKLEB.csv"
  )) %>%
    mutate(
      Taxon = gsub("#", "_", Taxon),
      Cluster = paste0("K", Cluster)
    ),
  read_csv(here(
    "data-raw/poppunk/strain_db_clustersESCO.csv"
  )) %>%
    mutate(
      Taxon = gsub("#", "_", Taxon),
      Cluster = paste0("E", Cluster)
    )
) -> btESBL_popPUNK

left_join(
  btESBL_sequence_sample_metadata,
  btESBL_popPUNK,
  by = c("lane" = "Taxon")
) -> btESBL_sequence_sample_metadata

use_data(btESBL_sequence_sample_metadata, overwrite = TRUE)

# contig clusters --------------------------------


#make summary df

filez <- list.files(here("data-raw/contig_clusters/"))

out <- list()
for (i in 1:length(filez)) {
  dftemp <-
    read_tsv(paste0(
      here("data-raw/contig_clusters/", filez[i])
    ))
  dftemp$gene <-
    str_extract(filez[i], "(?<=contigs-).*(?=\\.clust\\.)")
  out[[i]] <- dftemp
}

do.call(rbind, out) -> clusters

clusters %>%
  mutate(clstr_name = paste0(gene, ".", clstr),
         lane = gsub("^\\.", "", id),
         lane = gsub("\\..*$", "", lane),
         clstr_iden = as.numeric(gsub("%", "", clstr_iden)),
         clstr_cov = as.numeric(gsub("%", "", clstr_cov))) %>%
  select(-clstr) %>%
  mutate(species = case_when(
    lane %in% btESBL_coregene_tree_kleb$tip.label ~ "K. pneumoniae",
    lane %in% btESBL_coregene_tree_esco$tip.label ~ "E. coli",
    TRUE ~ NA_character_)
    )-> btESBL_contigclusters

use_data(btESBL_contigclusters, overwrite = TRUE)


# snpdist matrices ----------------------------------------------------

snpdists.e <- read_csv(here("data-raw/snpdists/esco_clean.full.filtered_polymorphic_sites-snp-dists.csv"))



names(snpdists.e)[1] <- "sample"
snpdists.e %>%
  filter(sample != "Reference") %>%
  select(-Reference) %>%
  rename_with(~ gsub("#", "_", .x)) %>%
  mutate(sample =  gsub("#", "_", sample)) ->
  btESBL_snpdists_esco

# kleb
snpdists.k <- read_csv(here("data-raw/snpdists/clean.full.filtered_polymorphic_sites_kleb_snp-dists.csv"))

names(snpdists.k)[1] <- "sample"
snpdists.k %>%
  filter(sample != "Reference") %>%
  select(-Reference) %>%
  rename_with(~ gsub("#", "_", .x)) %>%
  mutate(sample =  gsub("#", "_", sample)) ->
  btESBL_snpdists_kleb

use_data(btESBL_snpdists_esco, overwrite = TRUE)
use_data(btESBL_snpdists_kleb, overwrite = TRUE)

# ----------------------------------------------------------------




