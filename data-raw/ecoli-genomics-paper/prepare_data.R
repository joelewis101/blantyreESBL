# Prepare data for R package

#  accession numbers ---------------------------

library(blantyreESBL)
library(tidyverse)
library(devtools)

# load virulence determinents

vf <- read_csv(here("data-raw/all_dassim_esco_vf_ariba.csv"))


vf %>% 
  mutate(name = gsub("\\./", "", name),
         name = gsub("/report.tsv", "", name)) %>% 
  pivot_longer(-name, 
               names_to = c( "cluster", ".value"), 
               names_sep = "\\.") %>%
  filter(match == "yes") %>% 
  select(name, ref_seq) ->
  dassimEcoli_BTEcoli.virulence
  
use_data(dassimEcoli_BTEcoli.virulence, overwrite = TRUE)






btESBL_sequence_sample_metadata %>%
  filter(lane %in% btESBL_coregene_tree_esco$tip.label) %>%
  left_join(
    rbind(read_tsv(#"~/Documents/PhD/Thesis/bookdown/chapter_7/
      # checkm_quast/D1/transposed_report.tsv"
      here(
        "data-raw/QUAST_report1.tsv"
      )) ,
      read_tsv(#"~/Documents/PhD/Thesis/bookdown/chapter_7/
        #checkm_quast/D220190318/transposed_report.tsv"
        here(
          "data-raw/QUAST_report2.tsv"
        )),
      read_tsv(# "~/Documents/PhD/Thesis/bookdown/chapter_7/
        #  checkm_quast/D220190503/transposed_report.tsv"
        here(
          "data-raw/QUAST_report3.tsv"
        ))) %>%
      transmute(
        lane = gsub("\\.contigs_spades", "",  Assembly),
        number_of_contigs = `# contigs`,
        N50 = N50
      ),
    by = "lane"
  ) %>%
  rename(date_of_collection = data_date) %>%
  select(accession,
         lane,
         supplier_name,
         pid,
         date_of_collection,
         number_of_contigs,
         N50) %>%
  # add in phylogroups
  left_join(# rownames(df_hclusts) <- df_hclusts$Taxon
    left_join(read_csv(#"~/Documents/PhD/Thesis/bookdown/chapter_7/
      #phylogroup_and_mlst/mlst.csv"
      here(
        "data-raw/mlst.csv"
      )),
      read_csv(#"~/Documents/PhD/Thesis/bookdown/chapter_7/
        #phylogroup_and_mlst/phylogroups.csv"
        here(
          "data-raw/phylogroups.csv"
        )),
      by = c("lane" = "Lane")) %>%
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
    by = "lane") %>% 
# add in pathotypes
  left_join(
    dassimEcoli_BTEcoli.virulence %>% 
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
      select(name, Pathotype) %>% 
      unique() %>% 
      #pivot_wider(id_cols = name, values_from = pathotype, names_from = pathotype,
      #            values_fn = length, values_fill = NA_integer_) %>% 
      mutate(name = gsub("#", "_", name)),
    by = c("lane" = "name")
  ) %>% 
  left_join(
  # add in popPUNK clusters from global collection comparison
  read_csv(here("data-raw/popPUNK_clusters.csv")) %>% 
  mutate(
    Taxon = gsub("#","_",
                 gsub("\\..*$","", Taxon))),
  by = c("lane" = "Taxon")
  ) -> 
  dassimEcoli_BTEcoli.accession

use_data(dassimEcoli_BTEcoli.accession, overwrite = TRUE)

# add QRDR snps --------------------------

read_tsv(
  here("data-raw/QRDR_ariba_output.tsv")
  )  %>%
  rename(sample = "26141_1#222") %>% 
  filter(ref_name != "ref_name") %>% 
  mutate(sample = gsub("#", "_", sample)) %>%
  filter(sample %in% dassimEcoli_BTEcoli.accession$lane) %>%
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
  transmute(gene = ref_name,
            variant = ref_ctg_change,
            sample = sample) ->
  dassimEcoli_BTEcoli.QRDRmuts


dassimEcoli_CARD.QRDRmuts <- 
  read_csv(here("data-raw/2021-09-06_CARD-QRDR-mutations.csv"))

use_data(dassimEcoli_BTEcoli.QRDRmuts, overwrite = TRUE)
use_data(dassimEcoli_CARD.QRDRmuts, overwrite = TRUE)

# beta lactamase data
# downloaded 6 sept 2021


dassimEcoli_NCBI.betalactamases <- 
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
  )

use_data(dassimEcoli_NCBI.betalactamases, overwrite = TRUE)

# global metadata table --------------------


read_csv("data-raw/F1_genome_metadata.csv") ->
  dassimEcoli_Horesh.metadata

use_data(dassimEcoli_Horesh.metadata, overwrite = TRUE)

# TODO - merge in accessions

read_csv(here("data-raw/musicha_Ecoli_MetaData.csv")) %>%
  mutate(Lane = gsub("#", "_", Lane)) %>%
  left_join(
    left_join(read_csv(here( "data-raw/musicha_mlst.csv" )),
      read_csv( here( "data-raw/musicha_pgroup.csv" )),
      by = c("lane" = "V1")) %>%
      select(lane, ST, phylogroup) %>% 
      left_join( 
        read_csv(here( "data-raw/popPUNK_clusters.csv" )) %>%
          mutate(Taxon = gsub(
            "#", "_",
            gsub("\\..*$", "", Taxon))),
        by = c("lane" = "Taxon")),
    by = c("Lane" = "lane") 
  ) %>% 
  left_join(
    read_csv("data-raw/musicha_accession.csv") %>% 
      transmute(Accession = `Lane accession`,
                Lane = gsub("#","_", `Lane name`)),
    by = "Lane") -> dassimEcoli_Musicha.metadata


use_data(dassimEcoli_Musicha.metadata, overwrite = TRUE)

# global tree -----------------------

ape::read.tree(here("data-raw/IQTREE_globaltree.treefile")) -> 
  dassimEcoli_globaltree
phytools::midpoint.root(dassimEcoli_globaltree) -> dassimEcoli_globaltree

use_data(dassimEcoli_globaltree, overwrite = TRUE)
