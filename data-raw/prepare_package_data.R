# Load raw data and save it as .rda files for blantyreESBL package

# the load_phd_data scripts are from the Thesis repo https://github.com/joelewis101/thesis
#library(tidyverse)


source("/Users/joelewis/Documents/PhD/Thesis/bookdown/final_cleaning_scripts/load_PhD_data.R")

source("/Users/joelewis/Documents/PhD/Thesis/bookdown/final_cleaning_scripts/load_and_clean_lims.R")

library(tidyverse)
library(lubridate)
library(phytools)
library(devtools)
library(PopGenome)
library(here)
library(ape)

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

# longitudinal exposure (covariate) data ---------------------------------

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

# stan_data -------------------------------

readRDS("data-raw/stan_data_m2.rds") ->
  btESBL_stanmodeldata

use_data(btESBL_stanmodeldata , overwrite = TRUE)

read_csv("/Users/joelewis/Documents/PhD/Thesis/bookdown/data/stan_df.csv") %>%
  select(-`...1`) ->
  btESBL_modeldata

use_data(btESBL_modeldata, overwrite = TRUE)

readRDS("/Users/joelewis/Documents/PhD/Thesis/bookdown/chapter_9/stan_models/model_2/stanfit_m2.rds") ->
  btESBL_model2posterior

use_data(btESBL_model2posterior, overwrite = TRUE)

readRDS("/Users/joelewis/Documents/PhD/Thesis/bookdown/chapter_9/stan_models/model_1/stanfit_m1.rds") ->
  btESBL_model1posterior

use_data(btESBL_model1posterior, overwrite = TRUE)

# model3 -

list.files(here("data-raw/review_comment_work/models/"),
           pattern = "csv") -> csvfiles

rstan::read_stan_csv(
  paste0(here("data-raw/review_comment_work/models/"),
              csvfiles)
) -> btESBL_model3posterior


use_data(btESBL_model3posterior, overwrite = TRUE)

# simulations from posterior

btESBL_model2simulations <-
  readRDS(here("data-raw/btESBL_model2simulations.rda"))


use_data(btESBL_model2simulations, overwrite = TRUE)



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
    genus = "E. coli"),

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
    genus = "E. coli"),

klebs %>%
  mutate(name = gsub("\\./","", name),
         name = gsub("/report\\.tsv","",name)) %>%
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
    genus = "K. pneumoniae complex"
    )
) %>%
  mutate(name = gsub("\\./","", name),
         name = gsub("/report\\.tsv","",name),
         name = gsub("_1\\.fastq\\.gz", "", name)
  )  %>%
  rename("lane" = "name")  %>%
  mutate(lane = gsub("#", "_", lane)) %>%
  filter(lane %in% btESBL_sequence_sample_metadata$lane) %>%
  mutate(ref_seq = if_else(ref_seq == "TEM_95",
                        "TEM_1",
                        ref_seq)) ->
  btESBL_amrgenes

use_data(btESBL_amrgenes, overwrite = TRUE)

# QRDR resistance determinants


bind_rows(
  read_tsv(
    here("data-raw/ecoli-genomics-paper/QRDR_ariba_output.tsv")
  )  %>%
    rename(sample = "26141_1#222") %>%
    filter(ref_name != "ref_name") %>%
    mutate(sample = gsub("#", "_", sample)) %>%
    filter(sample %in% btESBL_sequence_sample_metadata$lane) %>%
    mutate(
      codon_posn = str_extract(ref_ctg_change, "[0-9]+"),
      codon_posn = as.numeric(codon_posn)
    ) %>%
    filter(ref_ctg_effect == "NONSYN") %>%
    filter(
      (ref_name == "GyrA" &
         codon_posn >= 67 & codon_posn <= 106) |
        (ref_name == "GyrB" &
           codon_posn >= 426 & codon_posn <= 464) |
        (ref_name == "ParC" &
           codon_posn >= 56 & codon_posn <= 108) |
        (ref_name == "ParE" &
           codon_posn >= 365 & codon_posn <= 525)
    ) %>%
    mutate(ref_name =  gsub("(^.{1})", '\\L\\1',
                            ref_name,
                            perl = TRUE)) %>%
    transmute(
      gene = ref_name,
      variant = ref_ctg_change,
      lane = sample,
      genus = "E. coli"
    ),
  # kleb
  read_tsv(
    here(
      "data-raw/kleb-genomics-paper/DASIM3_ariba_qrdr_snpcalls.tsv"
    )
  ) %>%
    rename(sample  = `34154_7#184`) %>%
    filter(ref_name != "ref_name") %>%
    mutate(sample = gsub("#", "_", sample)) %>%
    filter(sample %in% btESBL_sequence_sample_metadata$lane) %>%
    mutate(
      codon_posn = str_extract(ref_ctg_change, "[0-9]+"),
      codon_posn = as.numeric(codon_posn)
    ) %>%
    filter(ref_ctg_effect == "NONSYN") %>%
    filter(
      (ref_name == "GyrA" &
         codon_posn >= 67 & codon_posn <= 106) |
        (ref_name == "GyrB" &
           codon_posn >= 426 & codon_posn <= 464) |
        (ref_name == "ParC" &
           codon_posn >= 56 & codon_posn <= 108) |
        (ref_name == "ParE" &
           codon_posn >= 365 & codon_posn <= 525)
    ) %>%
    mutate(ref_name =  gsub("(^.{1})", '\\L\\1',
                            ref_name,
                            perl = TRUE)) %>%
    transmute(
      gene = ref_name,
      variant = ref_ctg_change,
      lane = sample,
      genus = "K. pneumoniae complex"
    )
) ->
  btESBL_qrdr_mutations

use_data(btESBL_qrdr_mutations, overwrite = TRUE)

# CARD described QRDR mutations -----------------------

btESBL_CARD_qrdr_mutations <-
  read_csv(
    here("data-raw/ecoli-genomics-paper/2021-09-06_CARD-QRDR-mutations.csv"))

use_data(btESBL_CARD_qrdr_mutations, overwrite = TRUE)

# NCBI beta-lactamase definition

btESBL_NCBI_phenotypic_bl <-
  read_tsv("https://ftp.ncbi.nlm.nih.gov/pathogen/betalactamases/Allele.tab") %>%
  rename_with( ~ tolower(gsub(" ", "_", .x))) %>%
  rename_with( ~ gsub("#", "", .x)) %>%
  mutate(
    allele_name = gsub("-", "_", allele_name),
    class = case_when(
      grepl("carbapenem-hydrolyzing",
            curated_gene_product_name) ~ "Carbapenemase",
      grepl("metallo-beta-lactamase",
            curated_gene_product_name) ~ "Carbapenemase",
      grepl(
        "extended-spectrum beta-lactamase",
        curated_gene_product_name
      ) ~ "ESBL",
      grepl("class C",
            curated_gene_product_name) ~ "AmpC",
      grepl("beta-lactamase",
            curated_gene_product_name) ~ "Penicillinase"
    )
  ) %>%
  select(allele_name, protein_accession_, nucleotide_accession_,
         gene_name, curated_gene_product_name, class)

use_data(btESBL_NCBI_phenotypic_bl, overwrite = TRUE)

# core gene trees --------------------------------------------------------

read.tree(here("data-raw/core_gene_trees/esco_core_gene_tree.treefile")) ->
  btESBL_coregene_tree_esco

read.tree(here("data-raw/core_gene_trees/kleb_core_gene_tree.treefile")) ->
  btESBL_coregene_tree_kleb

btESBL_coregene_tree_esco <- midpoint.root(btESBL_coregene_tree_esco)
btESBL_coregene_tree_kleb <- midpoint.root(btESBL_coregene_tree_kleb)

use_data(btESBL_coregene_tree_esco, overwrite = TRUE)
use_data(btESBL_coregene_tree_kleb, overwrite = TRUE)

# non ASC trees

read.tree(here(
  "data-raw/core_gene_trees/kleb_core_gene_tree_nonASC.treefile")) ->
  btESBL_coregene_tree_kleb_nonASC

read.tree(here(
  "data-raw/core_gene_trees/esco_core_gene_tree_nonASC.treefile")) ->
  btESBL_coregene_tree_esco_nonASC

btESBL_coregene_tree_kleb_nonASC <-
  midpoint.root(btESBL_coregene_tree_kleb_nonASC)

btESBL_coregene_tree_esco_nonASC <-
  midpoint.root(btESBL_coregene_tree_esco_nonASC)

#use_data(btESBL_coregene_tree_kleb_nonASC, overwrite = TRUE)

# global kleb tree -----------------------------------------------

btESBL_kleb_globaltree <-
  read.tree(
    here("data-raw/GLOBAL_core_gene_alignment_snp_sites.fa.treefile"))
midpoint.root(btESBL_kleb_globaltree) -> btESBL_kleb_globaltree

# global e coli tree

btESBL_ecoli_globaltree_noASC <-
  read.tree(
    here("data-raw/ecoli-genomics-paper/horesh_full_diversity/core_gene_alignment_snpsites.fasta.treefile"
         ))

midpoint.root(btESBL_ecoli_globaltree_noASC) ->
  btESBL_ecoli_globaltree_noASC

use_data(btESBL_ecoli_globaltree_noASC, overwrite = TRUE)

# use_data(btESBL_kleb_globaltree, overwrite = TRUE)

# non ASC treee

btESBL_kleb_globaltree_noASC <-
  read.tree(
    here("data-raw/kleb-genomics-paper/GLOBAL_core_gene_alignment_snp_sites_noASC.fa.treefile"))
midpoint.root(btESBL_kleb_globaltree_noASC) ->
  btESBL_kleb_globaltree_noASC


# malawi kleb tree -----------------------------

btESBL_kleb_malawi_allisolate_core_gene_tree <-
  read.tree(
    here("data-raw/MALAWI_core_gene_alignment_snp_sites.fa.treefile")
  )
midpoint.root(btESBL_kleb_malawi_allisolate_core_gene_tree) ->
  btESBL_kleb_malawi_allisolate_core_gene_tree
use_data(btESBL_kleb_malawi_allisolate_core_gene_tree, overwrite = TRUE)

# no ASC

btESBL_kleb_malawi_allisolate_core_gene_tree_noASC <-
  read.tree(
    here("data-raw/kleb-genomics-paper/MALAWI_core_gene_alignment_snp_sites_noASC.fa.treefile")
  )
midpoint.root(btESBL_kleb_malawi_allisolate_core_gene_tree_noASC) ->
  btESBL_kleb_malawi_allisolate_core_gene_tree_noASC
# use_data(btESBL_kleb_malawi_allisolate_core_gene_tree_noASC, overwrite = TRUE)

# overwrite ASC trees withnon ASC - kleb
btESBL_coregene_tree_kleb <- btESBL_coregene_tree_kleb_nonASC

use_data(btESBL_coregene_tree_kleb, overwrite = TRUE)

btESBL_kleb_globaltree <- btESBL_kleb_globaltree_noASC
use_data(btESBL_kleb_globaltree, overwrite = TRUE)

btESBL_kleb_malawi_allisolate_core_gene_tree <-
  btESBL_kleb_malawi_allisolate_core_gene_tree_noASC
use_data(btESBL_kleb_malawi_allisolate_core_gene_tree, overwrite = TRUE)

# and esco

btESBL_coregene_tree_esco <- btESBL_coregene_tree_esco_nonASC
use_data(btESBL_coregene_tree_esco, overwrite = TRUE)

btESBL_ecoli_globaltree_noASC -> btESBL_ecoli_globaltree

use_data(btESBL_ecoli_globaltree, overwrite = TRUE)

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
  sanger_metadata_file = "/Users/joelewis/Documents/Sanger/DASSIM3/data_raw/metadata/all_DASSIM_kleb_and_esco_metadata.csv",
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
  read_csv("data-raw/all_dassim_kleb_and_esco_accession.csv") %>%
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
  relocate(accession, .before = everything()) ->
  btESBL_sequence_sample_metadata

#  E coli virulence determinants --------------------------------------


vf <- read_csv(
  here("data-raw/ecoli-genomics-paper/all_dassim_esco_vf_ariba.csv"))


vf %>%
  mutate(name = gsub("\\./", "", name),
         name = gsub("/report.tsv", "", name)) %>%
  pivot_longer(-name,
               names_to = c( "cluster", ".value"),
               names_sep = "\\.") %>%
  filter(match == "yes") %>%
  select(name, ref_seq) ->
  btESBL_ecoli_virulence

use_data(btESBL_ecoli_virulence, overwrite = TRUE)


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
    "data-raw/poppunk/strain_db_clusters.csv"
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

# add assembly statistics --------------------------
btESBL_sequence_sample_metadata %>%
  left_join(
# e coli
rbind(read_tsv(#"~/Documents/PhD/Thesis/bookdown/chapter_7/
  # checkm_quast/D1/transposed_report.tsv"
  here(
    "data-raw/ecoli-genomics-paper/QUAST_report1.tsv"
  )) ,
  read_tsv(#"~/Documents/PhD/Thesis/bookdown/chapter_7/
    #checkm_quast/D220190318/transposed_report.tsv"
    here(
      "data-raw/ecoli-genomics-paper/QUAST_report2.tsv"
    )),
  read_tsv(# "~/Documents/PhD/Thesis/bookdown/chapter_7/
    #  checkm_quast/D220190503/transposed_report.tsv"
    here(
      "data-raw/ecoli-genomics-paper/QUAST_report3.tsv"
    ))) %>%
  transmute(
    lane = gsub("\\.contigs_spades", "",  Assembly),
    number_of_contigs = `# contigs`,
    N50 = N50
  ) %>%
  bind_rows(
    read_tsv(
      here(
        "data-raw/kleb-genomics-paper/DASSIM3_QUAST_transposed_report.tsv")
      ) %>%
      transmute(
        lane = gsub("\\.contigs_spades", "",  Assembly),
        number_of_contigs = `# contigs`,
        N50 = N50
      )
    ),
by = "lane"
) %>%
  # add MLST calls -----------------------------------------------
  #- e coli
  left_join(
    bind_rows(
    read_csv(here("data-raw/ecoli-genomics-paper/mlst.csv")) %>%
      mutate(ST = case_when(
        lane %in% c(
          # update those that were novel but now not
          "28099_1_144",
          "28099_1_18",
          "28099_1_159",
          "28099_1_249",
          "28099_1_23",
          "28099_1_249",
          "28099_1_330",
          "26141_1_141"
        ) ~ "9847",
        TRUE ~ ST
      )),
    # klebs
    read_tsv(
      here("data-raw/kleb-genomics-paper/DASSIM3_ariba_mlst_summary.tsv")
      ) %>%
      filter(ST != "ST") %>%
      mutate(ST = gsub("\\*", "", ST),
             lane = gsub("#","_", lane)) %>%
      select(ST, lane)
    ),
    by = "lane"
  ) %>%
  # add e coli phylogroup ------------------------------------------
  left_join(
    read_csv(here("data-raw/ecoli-genomics-paper/phylogroups.csv")) %>%
      transmute(lane = Lane,
                ecoli_phylogroup = Phylogroup),
    by = "lane") %>%
  # e coli pathotype
  left_join(
    btESBL_ecoli_virulence %>%
      filter(grepl("stx|ltcA|sta|eae|aatA|aggR|aaiC|ipaH|ipaD", ref_seq)) %>%
      mutate(ref_seq = gsub("_[0-9|A-Z|a-z]*$","", ref_seq)) %>%
      mutate(Pathotype = case_when(
        grepl("stx", ref_seq) & grepl("eae", ref_seq) ~ "EHEC",
        grepl("stx", ref_seq) ~ "STEC",
        grepl("eae", ref_seq) ~ "aEPEC/EPEC",
        grepl("aatA|aggR|aaiC", ref_seq) ~ "EAEC",
        grepl("ltcA|sta", ref_seq) ~ "ETEC",
        grepl("ipaH|ipaD", ref_seq) ~ "EIEC",
        TRUE ~ NA_character_ )) %>%
      transmute(lane = name,
                ecoli_pathotype = Pathotype) %>%
      unique() %>%
      mutate(lane = gsub("#", "_", lane)),
    by = "lane"
  ) %>%
  # add kleb k, o locus, virulence
  left_join(
    read_tsv(here("data-raw/kleb-genomics-paper/DASSIM3_kleborate.all.txt")) %>%
      mutate(strain = gsub("\\.contigs_spades","",strain),
             strain = gsub("#", "_", strain)) %>%
      filter(ST != "ST") %>%
      transmute(
        lane = strain,
        species = species,
        kleb_k_locus = K_locus,
        kleb_k_locus_confidence = K_locus_confidence ,
        kleb_o_locus = O_locus,
        kleb_o_locus_confidence = O_locus_confidence,
        kleb_YbST = YbST,
        kleb_CbST = CbST,
        kleb_AbST = AbST,
        kleb_SmST = SmST,
        kleb_rmpA = rmpA,
        kleb_rmpA2 = rmpA2),
    by = "lane")  %>%
  mutate(species =
           if_else(is.na(species),
                   "E. coli",
                   species)
  ) ->
  btESBL_sequence_sample_metadata


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
    ) -> btESBL_contigclusters

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


# plasmid replicons --------------------------------

esco.plasm1 <-
  read_csv(here("data-raw/plasmids/D1ariba_plasmidfinder_clustsmall.csv"))
esco.plasm2 <-
  read_csv(here("data-raw/plasmids/D2ariba_plasmidfinder_summary.csv"))


klebs.plasm <- read_csv("data-raw/plasmids/klebs-plasmidfinder-summary.csv")

dassim.klebs <- read_lines(here("data-raw/included_lanes/kleb_lanes_retained_following_qc.txt"))


bind_rows(

  esco.plasm1 %>%
    mutate(name = gsub("\\./","", name),
           name = gsub("_1\\..+$","", name),
           name = gsub("#","_", name)) %>%
    pivot_longer(-name,
                 names_to = c("cluster_name", ".value"),
                 names_sep = "\\.") %>%
    filter(match == "yes") %>%
    select(name, ref_seq) %>%
    mutate(species = "E. coli") %>%
    filter(name %in% btESBL_coregene_tree_esco$tip.label),
  esco.plasm2 %>%
    mutate(name = gsub("\\./","", name),
           name = gsub("/report.tsv","", name),
           name = gsub("#","_", name)) %>%
    pivot_longer(-name,
                 names_to = c("cluster_name", ".value"),
                 names_sep = "\\.") %>%
    filter(match == "yes") %>%
    select(name, ref_seq) %>%
    mutate(species = "E. coli") %>%
    filter(name %in% btESBL_coregene_tree_esco$tip.label),

  klebs.plasm %>%
    mutate(name = gsub("\\./","", name),
           name = gsub("/report\\.tsv","",name),
           name = gsub("#","_", name)) %>%
    filter(name %in% gsub("#","_",dassim.klebs)) %>%
    pivot_longer(-name,
                 names_to = c("cluster_name", ".value"),
                 names_sep = "\\.") %>%
    filter(match == "yes") %>%
    select(name, ref_seq) %>%
    mutate(species = "K. pneumoniae")
) %>%
  rename("lane" = "name") -> btESBL_plasmidreplicons


use_data(btESBL_plasmidreplicons, overwrite = TRUE)

# ST410 data ----------------------------

st410_metadata <-
  read_tsv(here("data-raw/ecoli-genomics-paper/st410/st410.tsv"))

st410_metadata %>%
  select(
    Uberstrain,
    Name,
    `Data Source(Accession No.;Sequencing Platform;Sequencing Library;Insert Size;Experiment;Status)`,
    `Source Niche`,
    `Source Details`,
    Country,
    `Collection Year`,
    ST
  ) %>%  separate_rows(
    `Data Source(Accession No.;Sequencing Platform;Sequencing Library;Insert Size;Experiment;Status)`,
    sep = ","
  ) %>%
  separate(
    `Data Source(Accession No.;Sequencing Platform;Sequencing Library;Insert Size;Experiment;Status)`,
    into = c(
      "accession",
      "platform",
      "library",
      "insert_size",
      "experiment"
    ),
    sep = ";"
  ) ->
  btESBL_ecoli_st410_metadata

use_data(btESBL_ecoli_st410_metadata, overwrite = TRUE)

# st410 plasmids -----------------

st410_plasm <-
  read_csv(
    here("data-raw/ecoli-genomics-paper/st410/st410_pf_ariba_summary.csv"))


st410_plasm %>%
  mutate(name = gsub("\\./", "", name),
         name = gsub("/report.tsv", "", name),
         name = gsub("_filtered","", name),
         name = gsub("#","_", name)) %>%
  pivot_longer(-name,
               names_to= c( "cluster", ".value"),
               names_sep = "\\.") %>%
  filter(match == "yes") %>%
  select(name, ref_seq) %>%
  mutate(ref_seq = gsub("\\..*$", "", ref_seq),
         ref_seq = gsub("_.*$","", ref_seq),
         ref_seq = gsub("^FIA", "IncFIA", ref_seq)) ->
  btESBL_ecoli_st410_plasmids

use_data(btESBL_ecoli_st410_plasmids, overwrite = TRUE)

# st410 tree -----------------------------------------------

read.tree(
  here(
    "data-raw/ecoli-genomics-paper/st410/clean_full.filtered_pollymorphic_sites.ref_removed.snpsites.fasta.treefile")) ->
  st410_tree

midpoint.root(st410_tree) -> btESBL_ecoli_globalst410_tree

btESBL_ecoli_globalst410_tree$tip.label <-
  gsub("_filtered", "", btESBL_ecoli_globalst410_tree$tip.label)

use_data(btESBL_ecoli_globalst410_tree, overwrite = TRUE)

# non-ASC tree

read.tree(
  here(
    "data-raw/ecoli-genomics-paper/st410/non_ASC_trees/clean.full.filtered_polymorphic_sites.fasta.treefile")) ->
  btESBL_ecoli_globalst410_tree_noASC

midpoint.root(btESBL_ecoli_globalst410_tree_noASC) -> btESBL_ecoli_globalst410_tree_noASC
btESBL_ecoli_globalst410_tree_noASC$tip.label <-
  gsub("_filtered","", btESBL_ecoli_globalst410_tree_noASC$tip.label)

use_data(btESBL_ecoli_globalst410_tree_noASC, overwrite = TRUE)

# st410 amr ---------------------------------------------------


amr.ariba410 <-
  read_csv(
    here(
      "data-raw/ecoli-genomics-paper/st410/st410_ariba_srst2_summary.csv"))

amr.ariba410 %>%
  mutate(name = gsub("\\./", "", name),
         name = gsub("/report.tsv", "", name),
         name = gsub("_filtered", "", name),
         name = gsub("#","_", name)) %>%
  pivot_longer(-name,
               names_to= c( "cluster", ".value"),
               names_sep = "\\.") %>%
  mutate(gene = sapply(str_split(ref_seq, "__"), function(x) x[3])) %>%
  filter(match == "yes") %>%
  mutate(gene = case_when(
    gene == "TEM_95" ~ "TEM_1",
    TRUE ~ gene
  )) %>%
  select(name, gene) ->
  btESBL_ecoli_st410_amr

use_data(btESBL_ecoli_st410_amr, overwrite = TRUE)

### st167


st167_metadata <-
  read_tsv(here("data-raw/ecoli-genomics-paper/st167/st167.tsv"))

st167_metadata %>%
  select(
    Uberstrain,
    Name,
    `Data Source(Accession No.;Sequencing Platform;Sequencing Library;Insert Size;Experiment;Status)`,
    `Source Niche`,
    `Source Details`,
    Country,
    `Collection Year`,
    ST
  ) %>%  separate_rows(
    `Data Source(Accession No.;Sequencing Platform;Sequencing Library;Insert Size;Experiment;Status)`,
    sep = ","
  ) %>%
  separate(
    `Data Source(Accession No.;Sequencing Platform;Sequencing Library;Insert Size;Experiment;Status)`,
    into = c(
      "accession",
      "platform",
      "library",
      "insert_size",
      "experiment"
    ),
    sep = ";"
  ) ->
  btESBL_ecoli_st167_metadata

use_data(btESBL_ecoli_st167_metadata, overwrite = TRUE)

# st167 plasmids -----------------

st167_plasm <-
  read_csv(
    here("data-raw/ecoli-genomics-paper/st167/st167_pf_ariba_summary.csv"))


st167_plasm %>%
  mutate(name = gsub("\\./", "", name),
         name = gsub("/report.tsv", "", name),
         name = gsub("_filtered","", name),
         name = gsub("#","_", name)) %>%
  pivot_longer(-name,
               names_to= c( "cluster", ".value"),
               names_sep = "\\.") %>%
  filter(match == "yes") %>%
  select(name, ref_seq) %>%
  mutate(ref_seq = gsub("\\..*$", "", ref_seq),
         ref_seq = gsub("_.*$","", ref_seq),
         ref_seq = gsub("^FIA", "IncFIA", ref_seq)) ->
  btESBL_ecoli_st167_plasmids

use_data(btESBL_ecoli_st167_plasmids, overwrite = TRUE)

# st167 tree -----------------------------------------------

read.tree(
  here(
    "data-raw/ecoli-genomics-paper/st167/clean.full.filtered_polymorphic_sites.ref_removed.snpsites.fasta.treefile")) ->
  st167_tree

st167_tree$tip.label <- gsub("_filtered","", st167_tree$tip.label)

midpoint.root(st167_tree) -> btESBL_ecoli_globalst167_tree

use_data(btESBL_ecoli_globalst167_tree, overwrite = TRUE)

# non ASC tree

read.tree(
  here("data-raw/ecoli-genomics-paper/st167/tree_no_ASC/clean.full.filtered_polymorphic_sites.fasta.treefile")
) -> btESBL_ecoli_globalst167_tree_noASC


midpoint.root(btESBL_ecoli_globalst167_tree_noASC) ->
  btESBL_ecoli_globalst167_tree_noASC

btESBL_ecoli_globalst167_tree_noASC$tip.label <-
  gsub("_filtered","", btESBL_ecoli_globalst167_tree_noASC$tip.label)

use_data(btESBL_ecoli_globalst167_tree_noASC, overwrite = TRUE)

# overwrite non ASC trees

btESBL_ecoli_globalst410_tree <- btESBL_ecoli_globalst410_tree_noASC

use_data(btESBL_ecoli_globalst410_tree, overwrite = TRUE)

btESBL_ecoli_globalst167_tree <- btESBL_ecoli_globalst167_tree_noASC

use_data(btESBL_ecoli_globalst167_tree, overwrite = TRUE)

# st167 amr ---------------------------------------------------


amr.ariba167 <-
  read_csv(
    here(
      "data-raw/ecoli-genomics-paper/st167/st167_ariba_srst2_summary.csv"))

amr.ariba167 %>%
  mutate(name = gsub("\\./", "", name),
         name = gsub("/report.tsv", "", name),
         name = gsub("_filtered", "", name),
         name = gsub("#","_", name)) %>%
  pivot_longer(-name,
               names_to = c( "cluster", ".value"),
               names_sep = "\\.") %>%
  mutate(gene = sapply(str_split(ref_seq, "__"), function(x) x[3])) %>%
  filter(match == "yes") %>%
  mutate(gene = case_when(
    gene == "TEM_95" ~ "TEM_1",
    TRUE ~ gene
  )) %>%
  select(name, gene) ->
  btESBL_ecoli_st167_amr

use_data(btESBL_ecoli_st167_amr, overwrite = TRUE)

# st131 amr and tree---------------------

read.tree(here("data-raw/ecoli-genomics-paper/st131/non_ASC_trees/gub_base.filtered_polymorphic_sites.fasta.treefile")) ->
  btESBL_ecoli_globalst131_tree
phytools::midpoint.root(btESBL_ecoli_globalst131_tree) -> btESBL_ecoli_globalst131_tree

btESBL_ecoli_globalst131_tree$tip.label <-
  gsub("_filtered", "", btESBL_ecoli_globalst131_tree$tip.label)

use_data(btESBL_ecoli_globalst131_tree, overwrite = TRUE)

bind_rows(
  read_csv(
    here("data-raw/ecoli-genomics-paper/st131/mbio.00644-19-st001.csv")
  ) %>%
    mutate(clade = case_when(
      `BAPS-1` %in% c(3) ~ "A",
      `BAPS-1` %in% c(2,4,5) ~ "B",
      `BAPS-1` %in% c(1) ~ "C",
      TRUE ~ "other")) %>%
    select(Accession_number, Year,Country, clade),
  btESBL_sequence_sample_metadata %>%
    transmute(Accession_number = lane,
              Year = year(data_date),
              Country = "Malawi") %>%
    filter(Accession_number %in% btESBL_ecoli_globalst131_tree$tip.label)
) %>%
  bind_rows(btESBL_ecoli_musicha_metadata %>%
              transmute(Accession_number = lane,
                        Year = Year,
                        Country = "Malawi")) %>%
    filter(Accession_number %in% btESBL_ecoli_globalst131_tree$tip.label) %>%
  as.data.frame() -> btESBL_ecoli_st131_metadata

use_data(btESBL_ecoli_st131_metadata, overwrite = TRUE)
# amr


read_csv(
  here(
    "data-raw/ecoli-genomics-paper/st131/non_ASC_trees/st131_ariba_srst2_summary.csv"
  )) %>%
  mutate(name = gsub("\\./", "", name),
         name = gsub("/report.tsv", "", name),
         name = gsub("_filtered", "", name),
         name = gsub("#","_", name)) %>%
  pivot_longer(-name,
               names_to = c( "cluster", ".value"),
               names_sep = "\\.") %>%
  mutate(gene = sapply(str_split(ref_seq, "__"), function(x) x[3])) %>%
  filter(match == "yes") %>%
  mutate(gene = case_when(
    gene == "TEM_95" ~ "TEM_1",
    TRUE ~ gene
  )) %>%
  select(name, gene) ->
  btESBL_ecoli_st131_amr

use_data(btESBL_ecoli_st131_amr, overwrite = TRUE)

# plasmids

st167_plasm <-
  read_csv(
    here("data-raw/ecoli-genomics-paper/st167/st167_pf_ariba_summary.csv"))


read_csv(
  here("data-raw/ecoli-genomics-paper/st131/non_ASC_trees/st131-ariba-plasm-summary.csv")
  ) %>%
  mutate(name = gsub("\\./", "", name),
         name = gsub("/report.tsv", "", name),
         name = gsub("_filtered","", name),
         name = gsub("#","_", name)) %>%
  pivot_longer(-name,
               names_to= c( "cluster", ".value"),
               names_sep = "\\.") %>%
  filter(match == "yes") %>%
  select(name, ref_seq) %>%
  mutate(ref_seq = gsub("\\..*$", "", ref_seq),
         ref_seq = gsub("_.*$","", ref_seq),
         ref_seq = gsub("^FIA", "IncFIA", ref_seq)) ->
  btESBL_ecoli_st131_plasmids

use_data(btESBL_ecoli_st131_plasmids, overwrite = TRUE)

# merge in clermonTyping phylogroups

ct <- read_tsv(
  here("data-raw/ecoli-genomics-paper/horesh_full_diversity/analysis_2021-11-03_110200_phylogroups.txt")
, col_names = FALSE) %>%
  select(X1,X5) %>%
  rename(lane = X1,
         ct_phylogroup = X5)

# to horesh

btESBL_ecoli_horesh_metadata %>%
  mutate(matcher  =
           gsub("\\.|fasta|fa|trimmed|spades|velvet|contigs_|_trimmed","",
                Assembly_name)) %>% pull(matcher)

left_join(
  btESBL_ecoli_horesh_metadata %>%
    mutate(matcher  =
             gsub("\\.|fasta|fa|trimmed|spades|velvet|contigs_|_trimmed","",
                  Assembly_name)),
  ct %>%
    mutate(lane =
             gsub("\\.|fasta|fa|trimmed|spades|velvet|contigs_|_trimmed","",
                  lane)),
  by = c("matcher" = "lane")) %>%
  mutate(Phylogroup =
           case_when(Phylogroup =="Not Determined" &
                       !is.na(ct_phylogroup) ~ ct_phylogroup,
                     TRUE ~ Phylogroup)) %>%
  select(-c(ct_phylogroup, matcher)) -> btESBL_ecoli_horesh_metadata

# dassim samples

btESBL_sequence_sample_metadata %>%
  left_join(
    ct %>%
      mutate(lane =
               gsub("\\.|fasta|fa|trimmed|spades|velvet|contigs_|_trimmed","",
                    lane),
             lane = gsub("#","_", lane))
  ) %>%
  mutate(ecoli_phylogroup = ct_phylogroup) %>%
  select(-ct_phylogroup) -> btESBL_sequence_sample_metadata

btESBL_ecoli_musicha_metadata %>%
  left_join(
    ct %>%
      mutate(lane =
               gsub("\\.|fasta|fa|trimmed|spades|velvet|contigs_|_trimmed","",
                    lane),
             lane = gsub("#","_", lane)) %>%
      mutate(phylogroup = ct_phylogroup) %>%
      select(-ct_phylogroup)
  ) -> btESBL_ecoli_musicha_metadata

use_data(btESBL_ecoli_horesh_metadata, overwrite = TRUE)

use_data(btESBL_sequence_sample_metadata, overwrite = TRUE)
use_data(btESBL_ecoli_musicha_metadata, overwrite = TRUE)


# load vfdb klebsiella output -------------------------------

vfdb <- read_csv("data-raw/kleb-genomics-paper/vfdb_core_ariba_summary.csv")

vfdb_genes <- read_lines("data-raw/kleb-genomics-paper/vfdb_core_gene_names.txt")

listout <- list()
for (i in 1:length(vfdb_genes)) {
  # print(i)
  strsplit(vfdb_genes[i], " |\\(|//)|>") -> split
  vfg <- grep("VFG", split[[1]], value = TRUE)
  vf <- grep("VF[0-9]", split[[1]], value = TRUE)
  desc <- str_extract(vfdb_genes[i], "(?<=\\) )[A-Za-z].*(?= \\[)")
  desc <- gsub("\\[.*$","", desc)
  if (length(vfg) == 0) {
    vfg <- NA_character_
  }
  if (length(vf) == 0) {
    vf <- NA_character_
  }
  listout[[i]] <-
    data.frame(
      vfg = vfg,
      vf = vf,
      description = desc
    )
}

do.call(rbind,listout) %>%
  mutate(vf = gsub("\\)\\]", "", vf)) %>%
  unique() ->
  vfdb_vf_lookup

vfdb %>%
  select(contains("name") |
           contains("match") |
           contains("assembled") |
           contains("ref_seq")) %>%
  pivot_longer(-name,
               names_to = c("cluster_name", ".value"),
               names_sep = "\\.") %>%
  filter(match == "yes") %>%
  mutate(ref_seq2 = sapply(strsplit(ref_seq, split = "\\."), function(x)
    x[1]),
    vfg = sapply(strsplit(ref_seq, split = "\\."), function(x)
      x[2])) %>%
  mutate(name = gsub("\\./|/report\\.tsv", "", name),
         name = gsub("#","_", name),
         vfg = gsub("_.*$","", vfg)) %>%
  select(name, ref_seq2, vfg) %>%
  mutate(vfg = if_else(is.na(vfg), ref_seq2, vfg)) -> vfdb_results


left_join(vfdb_results, vfdb_vf_lookup, by = "vfg")  -> vfdb_results

vfdb_results %>%
  mutate(ref_seq2 = gsub("afaE_I","afaE-I", ref_seq2)) %>%
  select(name,ref_seq2, vf, vfg) %>%
  # filter(
  #   vf != "VF0560",      #capsule
  #   vf != "VF0561",       #O type
  #   vf != "VF0564",    # yersinabactin
  #   vf != "VF0562" ,   # enterobactin
  #   vf != "VF0228", # enterobactin
  #   vf != "VF0564", # enterobactin
  #   !grepl("rmpA|rmpA2|ybt|iuc|iro|clb", ref_seq2)) %>%
  # select(name, vf, ref_seq2) %>%
  unique() %>%
  left_join(
    select(vfid_descriptions,
           VFID,
           VF_Name,
           VF_FullName),
    by = c("vf" = "VFID")) %>%
  rename(
    gene = ref_seq2) ->
  # mutate(VF_Name = if_else(
  #   !is.na(VF_FullName), VF_FullName,
  #   VF_Name)) %>%
  btESBL_kleb_malawi_vfdb_output

use_data(btESBL_kleb_malawi_vfdb_output, overwrite = TRUE)

# kleb global metadata ------------------

# load and clean global AMR - aim: ESBL vs not

amr.global <- read_csv(here("data-raw/kleb-genomics-paper/GLOBAL_ariba_amr.csv"))

amr.global %>%
  mutate(name = gsub("\\./", "", name),
         name = gsub("/report.tsv", "", name)) %>%
  pivot_longer(-name,
               names_to= c( "cluster", ".value"),
               names_sep = "\\.") %>%
  mutate(gene = sapply(str_split(ref_seq, "__"), function(x) x[3])) %>%
  filter(match == "yes") %>%
  mutate(gene = case_when(
    gene == "TEM_95" ~ "TEM_1",
    TRUE ~ gene
  )) %>%
  select(name, gene) %>%
  pivot_wider(names_from = gene,
              values_from = gene,
              values_fn = length,
              values_fill = 0) ->
  amr.global

amr.global %>%
  pivot_longer(-name, names_to = "gene") %>%
  filter(value == 1)  %>%
  left_join(select(btESBL_NCBI_phenotypic_bl, allele_name, class),
            by = c("gene" = "allele_name")) %>%
  mutate(class = case_when(
    gene == "SHV_12" ~ "ESBL",
    TRUE ~ class)
  ) %>%
  group_by(name) %>%
  summarise(ESBL =
              case_when(
                any(class == "ESBL") ~ "ESBL",
                TRUE ~ "0"
              )
  ) %>%
  # add accesssion
  left_join(
    bind_rows(
      read_csv(here("data-raw/kleb-genomics-paper/context_genomes/global_accession.csv")),
      read_csv(here("data-raw/kleb-genomics-paper/context_genomes/malawi_accession.csv"))
    ) %>%
      select(`Lane name`,
             `Sample accession`,
             `Lane accession`),
    by = c("name" = "Lane name")
  ) -> metadata_global



# merge in other data


musicha <- read_lines(here("data-raw/kleb-genomics-paper/context_genomes/musicha_klebs_list.txt"))
cornick <- read_lines(here("data-raw/kleb-genomics-paper/context_genomes/chathinka_kleb_lanes.txt"))
global <- read_lines(here("data-raw/kleb-genomics-paper/context_genomes/global_context_lanes.txt"))
kenya <- read_lines(here("data-raw/kleb-genomics-paper/context_genomes/kenya_lanes.txt"))


metadata_global %>%
  mutate(study =
           case_when(
             name %in% musicha ~ "musciha",
             name %in% cornick ~ "cornick",
             name %in% kenya ~ "kenya",
             name %in% global ~ "global",
             TRUE ~ "DASSIM"
           )) -> metadata_global

holt_metadata <-
  read_csv(here("data-raw/kleb-genomics-paper/context_genomes/holt_global_kleb_metadata.csv"))

holt_metadata %>%
  filter(!grepl("Kleb", File_ID)) %>%
  mutate(
    Year = if_else(
      Year_isolated == "Unknown",
      NA_real_,
      as.numeric(Year_isolated)),
    Infection_status =
           case_when(Source_Host == "Human" &
                       is.na(Infection_status) &
                       Sample_note %in% c(
                         "blood",
                         "sputum",
                         "urine",
                         "pus",
                         "bronchial alveolar lavage") ~ "Human_invasive",
                     Source_Host == "Human" &
                       is.na(Infection_status) &
                       Clinical_note == "Carriage" ~ "Human_carriage",
                     Source_Host == "Human" &
                       is.na(Infection_status) ~ "Unknown",
                     TRUE ~ Infection_status),
         sample_source = case_when(
           Source_Host == "Environmental" ~ "Environmental",
           Source_Host != "Human" ~ "Animal",
           TRUE ~ "Human"),
         isolate_type = case_when(
           Source_Host == "Human" &
             Infection_status == "Human_carriage" ~ "Carriage",
           Source_Host == "Human" &
             Infection_status %in% c(
               "Human_infection",
               "Human_invasive") ~ "Infection",
           Source_Host == "Human"  ~ NA_character_,
           TRUE ~ NA_character_)) -> holt_metadata

musicha_metadata <-
  read_csv(here("data-raw/kleb-genomics-paper/context_genomes/musciha_sample_metadata.csv"))

musicha_metadata %>%
  mutate(
    Year = Year,
    sample_source = "Human",
         isolate_type = case_when(
           Source == "RS" ~ "Carriage",
           TRUE~ "Infection")) -> musicha_metadata

left_join(
  metadata_global %>%
    mutate(name = gsub("#","_", name)),
  bind_rows(
    select(holt_metadata, File_ID,sample_source,isolate_type, Year),
    select(musicha_metadata, File_ID,sample_source,isolate_type, Year)
  ) %>%
    mutate(File_ID = gsub("#","_", File_ID)),
  by = c("name" = "File_ID")) %>%
  mutate(location =
           case_when(study %in% c(
             "musciha",
             "cornick",
             "DASSIM") ~ "Malawi",
             study == "kenya" ~ "Kenya",
             TRUE ~ "Global"),
         isolate_type =
           case_when(
             study %in% c("musicha","cornick", "kenya") ~ "Infection",
             study == "DASSIM" ~ "Carriage",
             TRUE ~ isolate_type),
         sample_source = case_when(
           study %in% c("musicha","cornick", "DASSIM","kenya") ~ "Human",
           TRUE ~ sample_source)
  ) %>%
  rename("Sample Source" = sample_source,
         "Isolate Type" = isolate_type) ->
  metadata_global

metadata_global %>%
  left_join(
    bind_rows(
      read_tsv(
        here("data-raw/kleb-genomics-paper/DASSIM3_kleborate.all.txt")) %>%
        mutate(strain = gsub("\\.contigs_spades","",strain)) %>%
        mutate(strain = gsub("\\.contigs_velvet", "", strain),
               ST = if_else(is.na(ST), "Unknown", ST)) %>%
        filter(ST != "ST") %>%
        select(
          strain,
          ST,
          YbST,
          CbST,
          AbST,
          SmST,
          rmpA,
          rmpA2,
          O_locus,
          K_locus,
          K_locus_confidence,
          O_locus_confidence
        ) %>%
        mutate(strain = gsub("#", "_", strain),
               ST = paste0("ST",ST)),
      read_tsv(
        here(
          "data-raw/kleb-genomics-paper/MALAWI_all_malawi_context_assemblies_kleborate_output.txt"
        )
      ) %>%
        mutate(strain = gsub("\\.contigs_velvet", "", strain),
               ST = if_else(is.na(ST), "Unknown", ST)) %>%
        filter(ST != "ST") %>%
        select(
          strain,
          ST,
          K_locus,
          K_locus_confidence,
          O_locus,
          O_locus_confidence,
          YbST,
          CbST,
          AbST,
          SmST,
          rmpA,
          rmpA2
        ) %>%
        mutate(ST = gsub("-.*$","", ST))
    ) %>%
      mutate(
        strain = gsub("#", "_", strain),
        ybt = if_else(YbST != "0", "1", "0"),
        clb = if_else(CbST != "0", "1", "0"),
        iuc = if_else(AbST != "0", "1", "0"),
        iro = if_else(SmST != "0", "1", "0"),
        rmpA = if_else(rmpA != "-", "1", "0"),
        rmpA2 = if_else(rmpA2 != "-", "1", "0")),
    by = c("name" = "strain")
  ) -> metadata_global

# add in malawi accessions
left_join(
  metadata_global,
  read_csv(
    here("data-raw/kleb-genomics-paper/DASSIM3_accessions_with_ERR.csv")
  ) %>%
    transmute(
      name = gsub("#", "_", `Lane name`),
      s_accession = `Sample name`,
      l_accession = `Lane accession`
    ),
  by = c("name" = "name")
) %>%
  left_join(
  btESBL_sequence_sample_metadata %>%
  transmute(
    lane = lane,
    year_dassim = lubridate::year(data_date),
    ST_ariba = if_else(ST == "Novel", "Novel",paste0("ST", ST))),
  by = c("name" = "lane")) %>%
  mutate(
    `Lane accession` =  if_else(is.na(`Lane accession`),
                 l_accession,
                 `Lane accession`),
    `Sample accession` =  if_else(is.na(`Sample accession`),
                 s_accession,
                 `Sample accession`),
    Year = if_else(
       study == "DASSIM",
       year_dassim,
       Year),
         ST = if_else(
           study == "DASSIM",
           ST_ariba,
           ST)) %>%
  select(-c(l_accession,s_accession, ST_ariba, -year_dassim)) ->
  metadata_global

metadata_global <- as.data.frame(metadata_global)
rownames(metadata_global) <- metadata_global$name
metadata_global$Malawi <- ifelse(metadata_global$location == "Malawi","1","0")
metadata_global$ST[metadata_global$ST == "0"] <- "Unknown"

metadata_global$ESBL <- as.character(metadata_global$ESBL)

btESBL_kleb_global_metadata <- metadata_global


use_data(btESBL_kleb_global_metadata, overwrite = TRUE)

# contig sens ax ----------------------------------

# function to load and prep data
source(here("data-raw/review_comment_work/contig_sens_ax/load_contig_sens_ax_data_fn.R"))
 btESBL_contigclusters_sensax <- load_contig_sens_ax_data()

 use_data(btESBL_contigclusters_sensax, overwrite = TRUE)

 # contig msa ----------------------------------------------------

 list.files(here("data-raw/review_comment_work/contig_msa"),
            pattern = "paf$",
            full.names = TRUE) -> paffiles
 list.files(here("data-raw/review_comment_work/contig_msa"),
            pattern = "paf$") -> paffile_names

 purrr::map(paffiles, read_tsv, col_select = 1:12,
     col_names =
       c("qname",
         "qlen",
         "qstart",
         "qend",
         "strand",
         "tname",
         "tlen",
         "tstart",
         "tend",
         "nmatch",
         "alen",
         "mapq")
 ) -> paf_file_list

 names(paf_file_list) <- gsub("\\.paf", "", paffile_names)

 btESBL_contigclusters_msa_paf_files <- paf_file_list

 list.dirs(here("data-raw/review_comment_work/contig_msa/alignments/"),
           recursive = FALSE) -> msa_dirs
 purrr::map(msa_dirs,
     PopGenome::readData,
     format = "fasta",
     include.unknown = TRUE) -> msa_list

 names(msa_list) <- gsub("\\.paf", "", paffile_names)

 btESBL_contigclusters_msa_alignments <- msa_list

use_data(btESBL_contigclusters_msa_paf_files, overwrite = TRUE)
use_data(btESBL_contigclusters_msa_alignments, overwrite = TRUE)


# BLAST results for msa ---------------------------------

blast_colnames <- c(
  "qseqid",
  "sseqid",
  "pident",
  "slen",
  "length",
  "mismatch",
  "gapopen",
  "qstart",
  "qend",
  "sstart",
  "send",
  "evalue",
  "bitscore"
)

bind_rows(
  read_csv(
    here(paste0(
      "data-raw/review_comment_work/plot_contigs/",
      "cluster_reps_srst2_blast.csv"
    )),
    col_names = blast_colnames
  ) %>%
    separate(sseqid,
             sep = "__",
             into = c(NA, "sseqid_group", "sseqid_gene", NA),
             remove = FALSE
    ) %>%
    mutate(type = "amr"),
  read_csv(
    here(paste0(
      "data-raw/review_comment_work/plot_contigs/",
      "cluster_reps_plasmidfinder_blast.csv"
    )),
    col_names = blast_colnames
  ) %>%
    mutate(sseqid = gsub("_.+$", "", sseqid)) %>%
    separate(sseqid,
             sep = "\\(",
             into = c("sseqid_group", "sseqid_gene"),
             remove = FALSE
    ) %>%
    mutate(sseqid_gene = gsub("\\)", "", sseqid_gene)) %>%
    mutate(type = "plasmid"),
  read_csv(
    here(paste0(
      "data-raw/review_comment_work/plot_contigs/",
      "cluster_reps_isfinder_blast.csv"
    )),
    col_names = blast_colnames
  ) %>%
    separate(sseqid,
             sep = "_",
             into = c("sseqid_gene", "sseqid_group", NA),
             remove = FALSE
    ) %>%
    mutate(type = "is")
) -> btESBL_contigclusters_msa_blastoutput

use_data(btESBL_contigclusters_msa_blastoutput, overwrite = TRUE)


# phenotypic AST data ------------------------------------------

btESBL_AST <-
  read_csv(
    "data-raw/review_comment_work/phenotypic_sens/ESBL_orgs.csv"
  )

btESBL_AST  %>%
  filter(profile_name == "DASSIM Culture") %>%
  transmute(
    supplier_name = sample_number,
    organism = organism,
    amikacin = `Amikacin 30`,
    chloramphenicol = `Chloramphenicol 30`,
    ciprofloxacin = `Ciprofloxacin 1`,
    cotrimoxazole = `Cotrimoxazole 25`,
    gentamicin = `Gentamicin 10`,
    meropenem = Meropenam) %>%
  unique() %>%
  filter(grepl("Escheric|ella pneum", organism)) %>%
  semi_join(
    btESBL_sequence_sample_metadata %>%
    mutate(
      organism = if_else(
        grepl("coli", species),
        "Escherichia coli",
        "Klebsiella pneumoniae")),
      by = c("supplier_name", "organism")
  )  %>%
  mutate(
    organism =
      case_when(
        organism == "Escherichia coli" ~ "E. coli",
        organism == "Klebsiella pneumoniae" ~ "KpSC")) %>%
  filter(!(is.na(amikacin) &
                 is.na(chloramphenicol) &
                 is.na(ciprofloxacin) &
                 is.na(cotrimoxazole) &
                 is.na(gentamicin) &
                 is.na(meropenem))) -> btESBL_AST

use_data(btESBL_AST, overwrite = TRUE)
