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
  filter(lane %in% btESBL_sequence_sample_metadata$lane) ->
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
    filter(name %in% dassim.klebs) %>%
    pivot_longer(-name,
                 names_to = c("cluster_name", ".value"),
                 names_sep = "\\.") %>%
    filter(match == "yes") %>%
    select(name, ref_seq) %>%
    mutate(species = "K. pneumoniae")
) %>%
  rename("lane" = "name") -> btESBL_plasmidreplicons


use_data(btESBL_plasmidreplicons, overwrite = TRUE)

### ST410 data ----------------------------

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
               names_to= c( "cluster", ".value"),
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
