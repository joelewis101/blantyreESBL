library(tidyverse)
library(ggtree)
library(here)
library(ape)
library(blantyreESBL)
library(lubridate)


read.tree(here("data-raw/ecoli-genomics-paper/st131/gub_base.filtered_polymorphic_sites.fasta.treefile")) ->
  btESBL_st131_tree
phytools::midpoint.root(btESBL_st131_tree) -> btESBL_st131_tree

df <- read_csv(here("data-raw/ecoli-genomics-paper/st131/mbio.00644-19-st001.csv")) %>%
  as.data.frame()
rownames(df) <- df$Accession_number
df$baps_1 <- as.character(df$`BAPS-1`)

btESBL_st131_tree$tip.label <- gsub("_filtered", "", btESBL_st131_tree$tip.label)

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
      filter(Accession_number %in% btESBL_st131_tree$tip.label)
    ) %>%
  as.data.frame() -> btESBL_st131_metadata

rownames(btESBL_st131_metadata) <- btESBL_st131_metadata$Accession_number

btESBL_st131_metadata %>%
  mutate(Malawi = if_else(Country == "Malawi",
                          "Malawi",
                          NA_character_))  -> btESBL_st131_metadata

(ggtree(btESBL_st131_tree) %<+% btESBL_st131_metadata +
  geom_tippoint(aes(color = Malawi)) +
  scale_color_manual(na.translate = FALSE, values = "red")) %>%
  gheatmap(select(btESBL_st131_metadata,clade))

(ggtree(treeio::tree_subset(btESBL_st131_tree, node= 1075, levels_back = 0)) %<+% btESBL_st131_metadata +
  geom_tippoint(aes(color = Malawi)) +
  scale_color_manual(na.translate = FALSE, values = "red")) %>%
  gheatmap(select(btESBL_st131_metadata,clade))  +
  geom_treescale()


