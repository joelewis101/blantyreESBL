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

  ggtree(btESBL_coregene_tree_esco ) |
    ggtree(btESBL_coregene_tree_esco_nonASC)

  ggtree(treeio::tree_subset(btESBL_kleb_globaltree_noASC, 734, levels_back = 0)) |
    ggtree(treeio::tree_subset(btESBL_coregene_tree_kleb_nonASC, 218, levels_back = 0))

  ggtree(
    treeio::tree_subset(btESBL_ecoli_globaltree_noASC, 3565, levels_back = 0),
    size = 0.3
  ) +
    geom_text(aes(label = node), hjust = -.3, size = 2)

  ggtree(
    btESBL_ecoli_globaltree_noASC,
    size = 0.3
  ) +
    geom_text(aes(label = node), hjust = -.3, size = 2)

  (
    (
      (
        (
          ggtree(btESBL_ecoli_globaltree_noASC, size = 0.4) %<+%
            (select(df_all, lane, ST) %>%
               mutate(
                 ST = if_else(ST == 410,
                              as.character("hilight"), NA_character_)
               ))  %>%
            gheatmap(
              select(df_all, `Year`),
              font.size = 4,
              width = 0.03,
              colnames_position = "top",
              color = NA,
              colnames_offset_y = 5,
              colnames_angle = 90,
              offset = 0.0006 + offset_add,
              hjust = 0
            )  +
            scale_fill_viridis_c(name = "Year", na.value = "white",
                                 guide = guide_colorbar(order = 2)) +
            new_scale_fill()
        ) %>%
          gheatmap(
            select(df_all, Continent),
            font.size = 4,
            width = 0.03,
            colnames_position = "top",
            color = NA,
            colnames_offset_y = 5,
            colnames_angle = 90,
            offset = 0.0012 + offset_add,
            hjust = 0
          ) +
          scale_fill_manual(values =
                              continent_cols,
                            name = "Continent",
                            guide = guide_legend(order = 3)) +
          # geom_tippoint(aes(color = ST), na.rm = TRUE, size = 0.5) +
          #   scale_color_manual(values  = "red", na.translate = FALSE) +
          new_scale_fill()
      ) %>%
        gheatmap(
          select(df_all, Country) %>%
            mutate(Country = if_else(
              Country == "Malawi",
              "Malawi", "Not Malawi"
            )),
          font.size = 4,
          width = 0.03,
          colnames_position = "top",
          color = NA,
          colnames_offset_y = 5,
          colnames_angle = 90,
          offset = 0.0018 + offset_add,
          hjust = 0
        ) +
        scale_fill_manual(values = country_cols,
                          name = "Country",
                          guide = guide_legend(order = 4)) +
        new_scale_fill()
    ) %>%
      gheatmap(
        select(df_all, Phylogroup) %>%
          mutate(Phylogroup = if_else(
            is.na(Phylogroup) ,"Unknown",
            Phylogroup)),
        width = 0.03,
        color = NA,
        font.size = 4,
        colnames_angle = 90,
        colnames_position = "top",
        colnames_offset_y = 3,
        hjust = 0,
        offset = 0 + offset_add
      ) +
      scale_fill_manual(values = pgroup_cols,
      name = "Phylogroup",
      na.translate = FALSE,
      guide = guide_legend(order = 1)) +
      new_scale_fill()
  ) +
    ylim(NA, 3900) +# + geom_tippoint(aes(color = ST)) +
    geom_treescale(y = 900, x = 0.005,offset = 10) +
    scale_color_manual(values = "red", na.translate = FALSE) -> p

  p + geom_hilight( node = 3550, alpha = 0.3,
                    fill = "firebrick",
                    color = "firebrick") +  # phylogroup A
    geom_cladelabel(node = 3550, label = "C", barsize = 0,offset = 0) +
    # node 1888 - st167
    # node 1848 - ST 617
    # node 1926 - ST 44
    # 5404
    geom_highlight(node = 5404, alpha = 0.3, fill = "firebrick",
                   color = "firebrick", extend = 0.0002) + # st 410 and context
    geom_cladelabel(node = 5404, label = "B", barsize = 0,offset = 0) +
    geom_highlight(node = 5652, alpha = 0.3, fill = "firebrick",
                   color = "firebrick") + #
    geom_cladelabel(node = 5652, label = "D", barsize = 0,offset = 0)
