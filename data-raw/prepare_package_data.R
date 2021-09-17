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


# add global context data

read_csv(here("data-raw/ecoli-genomics-paper/popPUNK_clusters.csv")) %>%
  transmute(
    lane = gsub("#","_",
                 gsub("\\..*$","", Taxon)),
    Cluster = Cluster) %>%
  semi_join(btESBL_sequence_sample_metadata,
            by = "lane") ->
  btESBL_ecoli_global_popPUNK_clusters

use_data(btESBL_ecoli_global_popPUNK_clusters,overwrite = TRUE)



read_csv("data-raw/ecoli-genomics-paper/F1_genome_metadata.csv") ->
  btESBL_ecoli_horesh_metadata

use_data(btESBL_ecoli_horesh_metadata, overwrite = TRUE)

# TODO - merge in accessions

read_csv(here("data-raw/ecoli-genomics-paper/musicha_Ecoli_MetaData.csv")) %>%
  mutate(Lane = gsub("#", "_", Lane)) %>%
  left_join(
    left_join(read_csv(here( "data-raw/ecoli-genomics-paper/musicha_mlst.csv" )),
              read_csv( here( "data-raw/ecoli-genomics-paper/musicha_pgroup.csv" )),
              by = c("lane" = "V1")) %>%
      select(lane, ST, phylogroup) %>%
      left_join(
        read_csv(here( "data-raw/ecoli-genomics-paper/popPUNK_clusters.csv" )) %>%
          mutate(Taxon = gsub(
            "#", "_",
            gsub("\\..*$", "", Taxon))),
        by = c("lane" = "Taxon")),
    by = c("Lane" = "lane")
  ) %>%
  left_join(
    read_csv("data-raw/ecoli-genomics-paper/musicha_accession.csv") %>%
      transmute(Accession = `Lane accession`,
                Lane = gsub("#","_", `Lane name`)),
    by = "Lane") %>%
  rename(lane = Lane) ->
  btESBL_ecoli_musicha_metadata


use_data(btESBL_ecoli_musicha_metadata, overwrite = TRUE)

# kleb global metadata


# load and clean global AMR - aim: ESBL vs not

amr.global <- read_csv(
  here("data-raw/kleb-genomics-paper/GLOBAL_ariba_amr.csv"))

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
  left_join(
    select(btESBL_NCBI_phenotypic_bl, allele_name, class),
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
      read_csv(
        here("data-raw/kleb-genomics-paper/context_genomes/global_accession.csv")),
      read_csv(
        here("data-raw/kleb-genomics-paper/context_genomes/malawi_accession.csv"))
    ) %>%
      select(`Lane name`,
             `Sample accession`),
    by = c("name" = "Lane name")
  )  -> metadata_global



# merge in other data


musicha <- read_lines(
  here("data-raw/kleb-genomics-paper/context_genomes/musicha_klebs_list.txt"))
cornick <- read_lines(
  here("data-raw/kleb-genomics-paper/context_genomes/chathinka_kleb_lanes.txt"))
global <- read_lines(
  here("data-raw/kleb-genomics-paper/context_genomes/global_context_lanes.txt"))
kenya <- read_lines(
  here("data-raw/kleb-genomics-paper/context_genomes/kenya_lanes.txt"))


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
  read_csv(
    here("data-raw/kleb-genomics-paper/context_genomes/holt_global_kleb_metadata.csv"))

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
  read_csv(here(
    "data-raw/kleb-genomics-paper/context_genomes/musciha_sample_metadata.csv"))

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
      btESBL_sequence_sample_metadata %>%
        transmute(
          strain = lane,
          ST = ST,
          YbST = kleb_YbST,
          CbST = kleb_CbST,
          AbST = kleb_AbST,
          SmST = kleb_SmST,
          rmpA = kleb_rmpA,
          rmpA2 = kleb_rmpA2,
          O_locus = kleb_o_locus,
          K_locus = kleb_k_locus,
          K_locus_confidence = kleb_k_locus_confidence,
          O_locus_confidence = kleb_o_locus_confidence
        ) %>%
        mutate(strain = gsub("#", "_", strain),
               ST = paste0("ST",ST)),
      read_tsv(
        here(
          "data-raw/kleb-genomics-paper/MALAWI_all_malawi_context_assemblies_kleborate_output.txt"
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

btESBL_kleb_global_metadata <- metadata_global


use_data(btESBL_kleb_global_metadata, overwrite = TRUE)

# global E. coli tree -----------------------

ape::read.tree(
  here("data-raw/ecoli-genomics-paper/IQTREE_globaltree.treefile")) ->
  btESBL_ecoli_globaltree

phytools::midpoint.root(btESBL_ecoli_globaltree) ->
  btESBL_ecoli_globaltree

use_data(btESBL_ecoli_globaltree, overwrite = TRUE)

# global klebtree -----------------



dassimKleb_trees.global <-
  read.tree(
    here("data-raw/kleb-genomics-paper/GLOBAL_core_gene_alignment_snp_sites.fa.treefile"))
midpoint.root(dassimKleb_trees.global) -> btESBL_kleb_globaltree

use_data(btESBL_kleb_globaltree, overwrite = TRUE)

#### malawi tree ---------------------------------------------

btESBL_kleb_malawi_allisolate_core_gene_tree <-
  read.tree(
    here("data-raw/kleb-genomics-paper/MALAWI_core_gene_alignment_snp_sites.fa.treefile")
  )
midpoint.root(btESBL_kleb_malawi_allisolate_core_gene_tree) ->
  btESBL_kleb_malawi_allisolate_core_gene_tree
use_data(btESBL_kleb_malawi_allisolate_core_gene_tree, overwrite = TRUE)
