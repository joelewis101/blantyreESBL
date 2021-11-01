library(tidyverse)
library(readxl)
library(here)
library(blantyreESBL)

# generate list of genomes to include from Gal's 10,000
# JL Oct 2021

# Strategy: include the top 500 she's selected. Include all isolates
# from popPUNK clusters size <=10, and sample remaining clusters down to
# 10 isolates

# load 500 clusters
metadata_top50 <- read_xlsx(
  here(
  "data-raw/ecoli-genomics-paper/st131/S1_table.xlsx"),
  skip = 4
)

metadata_others <- btESBL_ecoli_horesh_metadata %>%
  filter(PopPUNK > 51)

c(
  metadata_top50 %>%
    pull(Genome),
metadata_others %>%
  group_by(PopPUNK) %>%
  slice_sample(n=10) %>%
  mutate(Assembly_name = gsub("\\..*$","", Assembly_name)) %>%
  pull(Assembly_name)
) -> included_samples

write_lines(
  paste0(included_samples, ".gff"),
  here(
    paste0("data-raw/ecoli-genomics-paper/st131/",
  Sys.Date(),
  "horesh_genomes_sampled.txt"
  )))

# compare

gffs_filtered <- read_csv(
  here(
    "data-raw/ecoli-genomics-paper/horesh_full_diversity/gffs_filtered.txt"
  ),
  col_names = FALSE)



btESBL_sequence_sample_metadata %>%
  filter(species == "E. coli", ST == "131") -> st131

write_lines(
st131 %>%
  mutate(lane = gsub("_1_", "_1#", lane),
         lane = gsub("_2_", "_2#", lane)) %>%
  pull(lane),
here("data-raw/ecoli-genomics-paper/st131/dassim_st131.txt")
)

write_lines(
btESBL_ecoli_musicha_metadata %>%
  filter(ST == "131") %>%
  mutate(lane = gsub("_1_", "_1#", lane),
         lane = gsub("_8_", "_8#", lane)) %>%
  pull(lane),
here("data-raw/ecoli-genomics-paper/st131/musicha_st131.txt")
)

dep <- read_csv("data-raw/ecoli-genomics-paper/st131/depth_and_cov.csv")

dep$depth <- dep$depth_sum/dep$ref_bases
dep$cov <- dep$bases_mapped / dep$ref_bases

dep %>%
  select(lane, depth, cov) %>%
  pivot_longer((-lane)) %>%
  ggplot(aes(value)) +
  geom_histogram() +
  facet_wrap(~ name, scales = "free_x")

dep %>%
  ggplot(aes(depth, cov)) +
  geom_point() +
  geom_smooth()

  geom_histogram() +
  geom_vline(xintercept = 20)

dep %>%
  filter(depth <20) %>%
  pull(lane)
  group_by(depth < 20) %>%
  tally()

  write_lines(
  dep %>%
    filter(depth > 20) %>%
    mutate(lane = gsub("\\./|/snps.bam","",lane)
           ) %>%
    pull(lane),
  "data-raw/ecoli-genomics-paper/st131/st131_mappings_for_gubbins_following_qc.txt"
  )


  #### --------- reproduce analysis trees with non ASC trees

  ggtree(treeio::tree_subset(btESBL_coregene_tree_kleb, 210, levels_back = 0)) |
    ggtree(treeio::tree_subset(btESBL_coregene_tree_kleb_nonASC, 218, levels_back = 0))


  ggtree(btESBL_coregene_tree_kleb_nonASC) + geom_text(aes(label = node), hjust = -.3)

  ggtree(btESBL_kleb_globaltree_noASC) |
    ggtree(btESBL_kleb_globaltree)

  ggtree(treeio::tree_subset(btESBL_kleb_globaltree_noASC, 734, levels_back = 0)) |
    ggtree(treeio::tree_subset(btESBL_coregene_tree_kleb_nonASC, 218, levels_back = 0))

  ggtree(
    treeio::tree_subset(btESBL_kleb_globaltree_noASC, 734, levels_back = 0),
    size = 0.3
  ) +
    geom_text(aes(label = node), hjust = -.3, size = 2)

  btESBL_kleb_global_metadata %>%
    as.data.frame() ->
    df

  rownames(df) <- df$name

read_csv("data-raw/kleb-genomics-paper/context_genomes/holt_global_kleb_metadata.csv" ) %>%
  select(File_ID, ST) %>%
  mutate(ST = gsub("\\*|\\?","", ST)) %>%
  pivot_wider(id_cols = File_ID,
              names_from = ST,
              values_from = ST,
              values_fn = length)  %>%
  mutate(across(everything(), as.character)) %>%
  as.data.frame() -> df

rownames(df) <- df$File_ID

btESBL_kleb_global_metadata %>%
  select(name, ST) %>%
  filter(!is.na(ST)) %>%
  arrange(fct_infreq(ST)) %>%
  pivot_wider(id_cols = name,
              names_from = ST,
              values_from = ST,
              values_fn = length) %>%
  as.data.frame() -> mlst.onehot

rownames(mlst.onehot) <- mlst.onehot$name


ggtree( treeio::tree_subset(btESBL_kleb_malawi_allisolate_core_gene_tree_noASC,
                            372, levels_back = 0)) |
  ggtree(treeio::tree_subset(btESBL_kleb_malawi_allisolate_core_gene_tree, 372, levels_back = 0))

ggtree(
  # btESBL_kleb_malawi_allisolate_core_gene_tree_noASC
   treeio::tree_subset(btESBL_kleb_malawi_allisolate_core_gene_tree_noASC,
                       372, levels_back = 0),
   size = 0.3
) +
  geom_text(aes(label = node), hjust = -.3, size = 1)

ggtree(
  treeio::tree_subset(btESBL_kleb_malawi_allisolate_core_gene_tree_noASC, 372, levels_back = 0)
  ) %>%
  gheatmap(select(mlst.onehot, -name),
           width = 3,
           color = NA,
           font.size = 3,
           colnames_angle = 90,
           colnames_position = "top",
           colnames_offset_y = 0,
           hjust = 0,
           ) +
  ylim(NA, 350) -> p

# 627  = ST 218
(
  (
    ggtree(treeio::tree_subset(btESBL_kleb_malawi_allisolate_core_gene_tree_noASC, 372, levels_back = 0)) %>%
      gheatmap(
        select(dassimKleb_globalKleb.metadata, `Isolate Type`),
        width = 0.03,
        color = NA,
        font.size = 4,
        colnames_angle = 90,
        colnames_position = "top",
        colnames_offset_y = 5,
        hjust = 0,
        offset = 0.0002
      ) +
      scale_fill_manual(
        values =
          c(
            "Invasive" = brewer_pal(palette = "Set3")(6)[4],
            "Colonising" = brewer_pal(palette = "Set3")(6)[5]
          ),
        na.translate = FALSE,
        name = "Isolate\nType",
        guide = guide_legend(order = 2)
      ) +
      new_scale_fill()
  ) %>%   gheatmap(
    select(dassimKleb_globalKleb.metadata, ESBL)  %>%
      mutate(ESBL = if_else(ESBL == "ESBL", "Present", NA_character_)),
    width = 0.03,
    color = NA,
    font.size = 4,
    colnames_angle = 90,
    colnames_position = "top",
    colnames_offset_y = 5,
    hjust = 0,
    offset = 0
  )  +
    scale_fill_manual(
      values =
        c("Present" = brewer_pal(palette = "Set3")(6)[6]),
      na.translate = FALSE,
      name = "ESBL",
      guide = guide_legend(order = 1)
    ) +
    new_scale_fill()
) %>%   gheatmap(
  select(dassimKleb_globalKleb.metadata, ybt,clb,iuc,iro,rmpA,rmpA2) %>%
    mutate(across(everything(), ~
                    if_else(.x == "1", "Present", NA_character_))),
  width = 0.18,
  color = "lightgrey",
  font.size = 4,
  colnames_angle = 90,
  colnames_position = "top",
  colnames_offset_y = 5,
  hjust = 0,
  offset = 0.0005
)  +
  scale_fill_manual(values = c("Present" = "grey30"),
                    name = "Virulence\nGene", na.translate = FALSE,
                    guide = guide_legend(order = 3)) +
  ylim(NA, 370) +
  geom_cladelabel(node = 516, label = "ST268", align = TRUE, offset = - 0.0011) +
  geom_cladelabel(node = 629, label = "ST218", align = TRUE, offset = - 0.00163) +
  geom_cladelabel(node = 376, label = "ST14", align = TRUE, offset = - 0.00062) +
  geom_cladelabel(node = 412, label = "ST15", align = TRUE, offset = - 0.00074) +
  geom_cladelabel(node = 585, label = "ST340", align = TRUE, offset = - 0.00158) +
  geom_cladelabel(node = 523, label = "ST307", align = TRUE, offset = - 0.0013) +
  geom_treescale(x = 0.0003, y = 330, offset = 2, width = 0.001) -> malawi_tree_plot_final

malawi_tree_plot_final

