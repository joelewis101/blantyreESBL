library(tidyverse)
library(blantyreESBL)
library(readxl)
library(here)

library(ggtree)
library(ape)

st410_new <-
  read_xlsx(
    here("data-raw/ecoli-genomics-paper/st410-new/42003_2019_569_MOESM4_ESM.xlsx")
  )

st167_new_genomemed <-
  read_xls(
    here("data-raw/ecoli-genomics-paper/st167-new/13073_2019_699_MOESM2_ESM.xls")
  ) %>%
  filter(ST == 167)

st167_new_JAC <-
  read_xlsx(
    here("data-raw/ecoli-genomics-paper/st167-new/dky210_supplementary_tables.xlsx")
  ) %>%
  janitor::clean_names()

# almost all the JAC genomes are in the genomemed collection 

st167_new_genomemed %>%
  filter(!grepl("^G", Assembly) & !grepl("_", Assembly)) %>%
  pull(Assembly) %>%
  write_lines(
    here("data-raw/ecoli-genomics-paper/st167-new/st167_acc_new.txt")
    )

# dassim st167

btESBL_sequence_sample_metadata %>%
  filter(species == "E. coli", ST == 410) %>%
  nrow()

nrow(st167_new_genomemed)

st410_new %>%
  filter(Notes != "excluded" | is.na(Notes)) %>%
  pull(SRA)%>%
  write_lines(
    here("data-raw/ecoli-genomics-paper/st410-new/st410_acc_new.txt")
    )


st410_new %>%
  filter(Notes != "excluded" | is.na(Notes)) -> st410_new


#### st410 assemblies

st410_asm <- 
  read_csv(here("data-raw/ecoli-genomics-paper/st410-new/depth_and_cov.csv"))

st410_asm %>%
  mutate(dep = depth_sum/ bases_mapped,
         dep_whole_genome = depth_sum / ref_bases, 
         cov = bases_mapped / ref_bases) -> st410_asm 

ggplot(st410_asm, aes(fct_reorder(lane,dep), dep)) +
  geom_point() +
  geom_hline(yintercept = 20)

# Only 1 has sequencing depth below 20

st410_asm %>%
  filter(dep > 20) %>%
  mutate(lane = gsub("^\\.\\/|\\/snps.bam$","", lane)) %>%
  pull(lane) %>%
  write_lines(
    here("data-raw/ecoli-genomics-paper/st410-new/st410_assemblies_for_gubbins_following_qc.txt")
    )

  ### st167 assemblies


st167_asm <- 
  read_csv(here("data-raw/ecoli-genomics-paper/st167-new/depth_and_cov.csv"))

st167_asm %>%
  mutate(dep = depth_sum/ bases_mapped,
         dep_whole_genome = depth_sum / ref_bases, 
         cov = bases_mapped / ref_bases) -> st167_asm 

ggplot(st167_asm, aes(fct_reorder(lane,cov), cov)) +
  geom_point() 


sum(st167_asm$dep < 20)

# Only 2  has sequencing depth below 20

st167_asm %>%
  filter(dep > 20) %>%
  mutate(lane = gsub("^\\.\\/|\\/snps.bam$","", lane)) %>%
  pull(lane) %>%
  write_lines(
    here("data-raw/ecoli-genomics-paper/st167-new/st167_assemblies_for_gubbins_following_qc.txt")
    )



read.tree("data-raw/ecoli-genomics-paper/st167-new/clean.full.filtered_polymorphic_sites.fasta.treefile") -> st167tree

phytools::midpoint.root(st167tree) -> st167tree

data.frame(
           taxa = st167tree$tip.label,
country = if_else(grepl("^2", st167tree$tip.label ), "Malawi", NA_character_)) -> dd

ggtree(st167tree) %<+% dd + 
  geom_tippoint(aes(color = country))
