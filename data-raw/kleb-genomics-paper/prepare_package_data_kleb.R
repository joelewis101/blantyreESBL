# load and clean metadata and save

library(tidyverse)
library(here)
library(ape)
library(phytools)
library(devtools)

# Assembly stats, species, MLST, K and O type ---------------------------

# load data nnd merge
# This study Kleb diversity
left_join(
  read_csv(here("data-raw/DASSIM3_accession.csv")) %>%
    select("Sample accession","Lane name"),

  # add QUAST assembly stats
  read_tsv(here("data-raw/DASSIM3_QUAST_transposed_report.tsv")) %>%
    mutate(Assembly = gsub("\\.contigs_spades","", Assembly),
           Assembly = gsub("_7_","_7#", Assembly)) %>%
    select(Assembly, `# contigs`, N50, `Total length`),
  by = c("Lane name" = "Assembly")) %>%
  rename("No. contigs" = "# contigs")   %>%

  # add ARIBA mlst calls
  left_join(
    read_tsv(here("data-raw/DASSIM3_ariba_mlst_summary.tsv")) %>%
      filter(ST != "ST") %>%
      mutate(ST = gsub("\\*", "", ST)) %>%
      select(ST, lane),
    by = c("Lane name" = "lane")
  ) %>%

  #add kleborate O-, K- type, species, virulence determinents
  left_join(
    read_tsv(here("data-raw/DASSIM3_kleborate.all.txt")) %>%
      mutate(strain = gsub("\\.contigs_spades","",strain)) %>%
      filter(ST != "ST") %>%
      select(
        strain,
        species,
        species_match,
        K_locus,
        K_locus_confidence,
        O_locus,
        O_locus_confidence,
        YbST,CbST,AbST,SmST,rmpA,rmpA2),
    by = c("Lane name" = "strain")) ->
  dassimKleb_BTKleb.diversity

use_data(dassimKleb_BTKleb.diversity, overwrite = TRUE)


# Malawi tree -----------------------------------------------

dassimKleb_trees.BTcarriage <-
  read.tree(here("data-raw/DASSIM3core_gene_alignment_snp_sites.fa.treefile"))
midpoint.root(dassimKleb_trees.BTcarriage) -> dassimKleb_trees.BTcarriage

use_data(dassimKleb_trees.BTcarriage, overwrite = TRUE)

# AMR ----------------------------------------------------

amr.ariba <- read_csv(here("data-raw/DASSIM3-ariba-report.csv"))


# get beta lactamase phenotype and tidy

bls <- read_tsv(
  "https://ftp.ncbi.nlm.nih.gov/pathogen/betalactamases/Allele.tab"
)

bls %>%
  rename_with(~ tolower(gsub(" ", "_", .x))) %>%
  rename_with(~ gsub("#", "", .x)) %>%
  mutate(allele_name = gsub("-", "_", allele_name),
         class = case_when(
           grepl("carbapenem-hydrolyzing",
                 curated_gene_product_name) ~ "CPE",
           grepl("metallo-beta-lactamase",
                 curated_gene_product_name) ~ "CPE",
           grepl("extended-spectrum beta-lactamase",
                 curated_gene_product_name) ~ "ESBL",
           grepl("class C",
                 curated_gene_product_name) ~ "AmpC",
           grepl("beta-lactamase",
                 curated_gene_product_name) ~ "Penicillinase"
         )) ->
  bls


# get the ariba SNP for QRDR calls and parse

read_tsv(here("data-raw/DASIM3_ariba_qrdr_snpcalls.tsv")) -> amr.qrdr

names(amr.qrdr)[length(names(amr.qrdr))] <- "sample"
subset(amr.qrdr, ref_name != "ref_name") -> amr.qrdr

# rembember the weird and contam lanes are in this data - remove

include.lanes <-
  read_lines(
    here("data-raw/DASSIM3_lanes_retained_following_qc.txt"
  ))

# tidy up - restrict to QRDR (includes whole genes)

amr.qrdr %>%
  filter(sample %in% include.lanes) %>%
  mutate(codon_posn = str_extract(amr.qrdr$ref_ctg_change, "[0-9]+"),
         codon_posn = as.numeric(codon_posn)) %>%
  filter(ref_ctg_effect == "NONSYN") ->
  amr.qrdr

# refs https://www.frontiersin.org/articles/10.3389/fmicb.2015.01355/full for
# gyra and pac C
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3299416/ gyrB


amr.qrdr <- subset(
  amr.qrdr,
  (ref_name == "GyrA" &
     codon_posn >= 67 & codon_posn <= 106) |
    (ref_name == "GyrB" &
       codon_posn >= 426 & codon_posn <= 464) |
    (ref_name == "ParC" &
       codon_posn >= 56 & codon_posn <= 108) |
    (ref_name == "ParE" &
       codon_posn >= 365 & codon_posn <= 525)
)

# make one hot coded df

# tidy up other df , merge in QRDR
# remember TEM-95 is misnamed as TEM-1 - change

amr.ariba %>%
  mutate(name = gsub("\\./", "", name),
         name = gsub("/report.tsv", "", name)) %>%
  filter(name %in% include.lanes) %>%
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
  bind_rows(
  amr.qrdr %>%
  filter(grepl("87|82|S80", ref_ctg_change)) %>%
    # have checked that these are the only described ones in this df
  select(ref_name, sample) %>%
  unique() %>%
  dplyr::rename(
    gene = ref_name,
    name = sample
  )) %>%
  pivot_wider(names_from = gene,
              values_from = gene,
              values_fn = length,
              values_fill = 0) ->
  amr.ariba


# add class ro which resistnce is conferred, make long df

quinolone <- "Par|Gyr|Par|Qnr|Qep|Nor|GyrA|GyrB|ParC|ParE"
tetracycline<- "Tet"
sulphonamide <- "Sul"
aminoglycoside <- "Str|Aad|Aac|Aph|Rmt|APH"
streptothricin <- "Sat"
macrolide <- "Mph|Mdf|Erm|Ere"
fosfomycin <- "Fos"
chloramphenicol <- "Cat|FloR|Cml"
trimethoprim <- "Dfr"
rifampicin <- "Arr"
ESBL <- "SHV_12"
penicillinase <- "OKP|SCO|LEN|LAP|AmpH"
efflux <- "Oqx"


amr.ariba %>%
  pivot_longer(-name, names_to = "gene") %>%
#  filter(gene != "AmpH") %>%
  left_join(select(bls, allele_name, class),
            by = c("gene" = "allele_name")) %>%
  mutate(class = case_when(
    str_detect(gene, quinolone) ~ "Quinolone",
    str_detect(gene, tetracycline) ~ "Tetracycline",
    str_detect(gene, sulphonamide) ~ "Sulphonamide",
    str_detect(gene, aminoglycoside) ~ "Aminoglycoside",
    str_detect(gene, streptothricin) ~ "Streptothricin",
    str_detect(gene, macrolide) ~ "Macrolide",
    str_detect(gene, fosfomycin) ~ "Fosfomycin",
    str_detect(gene, chloramphenicol) ~ "Chloramphenicol",
    str_detect(gene, rifampicin) ~ "Rifampicin",
    str_detect(gene,trimethoprim) ~ "Trimethoprim",
    str_detect(gene,ESBL) ~ "ESBL",
    str_detect(gene,penicillinase) ~ "Penicillinase",
    str_detect(gene,efflux) ~ "Efflux",
    TRUE ~ class
  )) %>%
  filter(value == 1) %>%
  select(-value) ->
  dassimKleb_BTKleb.amr

use_data(dassimKleb_BTKleb.amr, overwrite = TRUE)

# Global metadata -------------------------------------------------------

# global kleb tree

btESBL_kleb_globaltree <-
  read.tree(
    here("data-raw/GLOBAL_core_gene_alignment_snp_sites.fa.treefile"))
midpoint.root(btESBL_kleb_globaltree) -> btESBL_kleb_globaltree

use_data(btESBL_kleb_globaltree, overwrite = TRUE)

# non ASC treee

btESBL_kleb_globaltree_noASC <-
  read.tree(
    here("data-raw/kleb-genomics-paper/GLOBAL_core_gene_alignment_snp_sites_noASC.fa.treefile"))
midpoint.root(btESBL_kleb_globaltree_noASC) ->
  btESBL_kleb_globaltree_noASC

use_data(btESBL_kleb_globaltree_noASC, overwrite = TRUE)

# load and clean global AMR - aim: ESBL vs not

amr.global <- read_csv(here("data-raw/GLOBAL_ariba_amr.csv"))

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
  left_join(select(bls, allele_name, class),
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
      read_csv(here("data-raw/context_genomes/global_accession.csv")),
      read_csv(here("data-raw/context_genomes/malawi_accession.csv"))
    ) %>%
      select(`Lane name`,
             `Sample accession`),
    by = c("name" = "Lane name")
  ) %>%
  relocate("Sample accession", before = everything()) -> metadata_global



# merge in other data


musicha <- read_lines(here("data-raw/context_genomes/musicha_klebs_list.txt"))
cornick <- read_lines(here("data-raw/context_genomes/chathinka_kleb_lanes.txt"))
global <- read_lines(here("data-raw/context_genomes/global_context_lanes.txt"))
kenya <- read_lines(here("data-raw/context_genomes/kenya_lanes.txt"))


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
  read_csv(here("data-raw/context_genomes/holt_global_kleb_metadata.csv"))

holt_metadata %>%
  filter(!grepl("Kleb", File_ID)) %>%
  mutate(Infection_status =
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
           TRUE ~ NA_character_)) ->holt_metadata

musicha_metadata <-
  read_csv(here("data-raw/context_genomes/musciha_sample_metadata.csv"))

musicha_metadata %>%
  mutate(sample_source = "Human",
         isolate_type = case_when(
           Source == "RS" ~ "Carriage",
           TRUE~ "Infection")) -> musicha_metadata

left_join(
  metadata_global %>%
    mutate(name = gsub("#","_", name)),
  bind_rows(
    select(holt_metadata, File_ID,sample_source,isolate_type),
    select(musicha_metadata, File_ID,sample_source,isolate_type)
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
    dassimKleb_BTKleb.diversity %>%
      rename(strain = "Lane name") %>%
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
        "data-raw/MALAWI_all_malawi_context_assemblies_kleborate_output.txt"
      )
    ) %>%
      mutate(strain = gsub("\\.contigs_velvet", "", strain)) %>%
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
      rmpA2 = if_else(rmpA2 != "-", "1", "0"))
    ,
    by = c("name" = "strain")
  ) -> metadata_global


metadata_global <- as.data.frame(metadata_global)
rownames(metadata_global) <- metadata_global$name
metadata_global$Malawi <- ifelse(metadata_global$location == "Malawi","1","0")

metadata_global$ESBL <- as.character(metadata_global$ESBL)

dassimKleb_globalKleb.metadata <- metadata_global


use_data(dassimKleb_globalKleb.metadata, overwrite = TRUE)

#### malawi tree ---------------------------------------------

dassimKleb_trees.malawi <-
  read.tree(
    here("data-raw/MALAWI_core_gene_alignment_snp_sites.fa.treefile")
  )
midpoint.root(dassimKleb_trees.malawi) -> dassimKleb_trees.malawi
use_data(dassimKleb_trees.malawi, overwrite = TRUE)

