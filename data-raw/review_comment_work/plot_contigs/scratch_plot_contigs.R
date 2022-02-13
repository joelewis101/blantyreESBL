# get contig rep clusters, make list

library(blantyreESBL)
library(tidyverse)
library(here)
library(IRanges)
library(gggenes)
library(ggrepel)
library(viridis)
library(ggnewscale)

btESBL_contigclusters %>%
  filter(clstr_rep == 1) %>%
  pull(id) -> cluster_reps

write_lines(
  cluster_reps,
  here("data-raw/review_comment_work/plot_contigs/cluster_reps.txt")
)

  # load all blast outputs into one df -------------

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
) -> blast_all


# functions for cleaning up blast output for plotting -----------------

# this functin will merge all overlapping matches
# then label them with range_group id
add_range_group_ids <- function(qstart, qend, sseqid) {
  ir <- IRanges(qstart, qend, names = sseqid)
  range_group <- subjectHits(findOverlaps(ir, reduce(ir)))
  return(range_group)
}



# take each range_group df
# if only one top bitscore
# select that_
# if ties - selct sseqid_group if they are sll one group
# if they are different groups, give a warning

pick_best_fitting_gene <- function(df) {
  # restrict only to top fittinmg bitscoire
  df %>%
    ungroup() %>%
    filter(bitscore == max(bitscore)) -> df
  if (nrow(df) == 1) {
    # if only one - set gene to sseqid gene
    df %>%
      mutate(
        gene = sseqid_gene,
        duplicate_flag = 0
      ) %>%
      select(-c(sseqid_group, sseqid_gene)) -> df
  } else {
    # if more than one top bitscore
    if (length(unique(df$sseqid_group)) == 1) {
      # if they are all same group, use group
      df %>%
        arrange(desc(pident), desc(length)) %>%
        dplyr::slice(n = 1) %>%
        mutate(
          gene = paste0(sseqid_group, "*"),
          duplicate_flag = 0
        ) %>%
        select(-c(sseqid_group, sseqid_gene)) -> df
    } else {
      df %>%
        group_by(sseqid_group) %>%
        arrange(sseqid_group, desc(pident), desc(length)) %>%
        dplyr::slice(n = 1) %>%
        mutate(
          gene = paste0(sseqid_group, "*"),
          duplicate_flag = 1
        ) %>%
        ungroup() %>%
        select(-c(sseqid_group, sseqid_gene)) -> df
      warning(
        paste0(
          "Watch out! qseqid: ", df$qseqid, " range_group: ",
          df$range_group, " have multiple genes in one location.\n"
        )
      )
    }
  }
  return(df)
}

# tidy up and plot -----------------------------------

blast_all %>%
  group_by(qseqid, type) %>%
  mutate(
    range_group = paste0(type,
                         add_range_group_ids( qstart, qend, sseqid ))
  ) -> blast_all

# test

blast_all %>%
  filter(qseqid == ".26141_1_284.27") -> test

pick_best_fitting_gene(filter(test, range_group == "is2"))

blast_all %>%
  ungroup() %>%
  group_by(qseqid, range_group) %>%
  filter(sseqid != "ISEc9_IS1380_unknown") %>%
  do(pick_best_fitting_gene(.)) -> blast_all_processed

# plot

blast_all_processed %>%
  filter(qseqid == ".28099_2_219.36") -> test

ggplot(test, aes(xmin = qstart, xmax = qend, y = qseqid,
                 fill = type,
                 label = gene)) +
  geom_gene_arrow(arrowhead_height = unit(7, "mm"),
                  arrowhead_width = unit(1, "mm"),
                  arrow_body_height = unit(7,"mm"),
                  alpha = 0.5) +
  # geom_blank(data = dummies1) +
  # facet_wrap(~ clust.id, scales = "free", ncol =1) + #theme_genes() +
  #theme(legend.position="none") +#  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(1, "mm")) +
  geom_text_repel(aes(x= qstart)) +
  theme_genes() +  theme(legend.position = "bottom") +
  scale_fill_discrete()

# try a plot added in to

#.28099_1_127.68 is CTXM15.123

paf_df <- btESBL_contigclusters_msa_paf_files$CTX_M_15.123
genes_df <- filter(blast_all_processed, qseqid == ".28099_1_127.68")


   paf_df %>%
    mutate(
      qstart_coord = if_else(
        strand == "+",
        tstart - qstart,
        tstart - (qlen - qend)
      ),
      qend_coord = qstart_coord + qlen,
      qcov_start_coord = if_else(
        strand == "+",
        qstart_coord + qstart,
        qstart_coord +
          (qlen - qend)
      ),
      qcov_end_coord = if_else(
        strand == "+",
        qstart_coord + qend,
        qstart_coord + (qlen - qstart)
      )
    ) %>%
    group_by(qname) %>%
    arrange(qlen, qname, alen) %>%
    mutate(
      map_id = row_number(),
      qname_map = paste0(qname, "_", map_id),
      map_type = if_else(alen == max(alen),
        "Primary",
        "Secondary"
      )
    ) %>%
    ungroup() %>%
    mutate(axis_scale_pos = if_else(
      map_type == "Primary", 0.8, 0.6
    )) %>%
    arrange(desc(qlen), desc(qname), desc(alen)) %>%
    mutate(axis_scale_pos = cumsum(axis_scale_pos)) -> df
  # fudgey axis positions - continuous scale that
  # will later be changed to 'discrete'
  df %>%
    mutate(axis_scale_pos = nrow(df) - axis_scale_pos) %>%
    group_by(qname) %>%
    mutate(axis_location = mean(axis_scale_pos)) -> df
# Add in reference
  bind_rows(
    data.frame(
      qname = "Cluster Reference",
      qstart_coord = 0,
      qend_coord = unique(df$tlen),
      qcov_start_coord = 0,
      qcov_end_coord = unique(df$tlen),
      axis_scale_pos = max(df$axis_scale_pos + 1),
      axis_location = max(df$axis_scale_pos + 1),
      map_type = "Primary",
      reference = "Reference"
    ),
    df
  ) %>%
    mutate(
      reference =
        case_when(
          is.na(reference) ~ strand,
          TRUE ~ reference
        )
    ) -> df

  genes_df %>%
    mutate(
      qname = "Cluster Reference",
      reference = NA
      ) -> genes_df


  # plot alignments

  df %>%
    mutate(
      across(starts_with("axis"),
                         ~ if_else(qname == "Cluster Reference",
                                   .x + 1, .x)
      )) -> df
  df %>%
    mutate(reference =
             if_else(strand %in% c("+", "-"),
                     "Covered", strand)) %>%
  ggplot(aes(
    xmin = qstart_coord,
    xmax = qend_coord,
    y = 1:length(qname_map),
    ymin = axis_scale_pos - 0.2,
    ymax = axis_scale_pos + 0.2,
    fill = reference,
    linetype = map_type
  )) +
    geom_rect(aes(
      xmin = qcov_start_coord,
      xmax = qcov_end_coord
    ),
    color = NA,
    size = 1,
    ) +
     geom_rect(
       color = "black",
       fill = "white",
       size = 0.5,
       alpha = 0
     ) +
    scale_y_continuous(
      labels = df$qname,
      breaks = df$axis_location
    ) +
    scale_x_continuous(labels = function(x) paste0(x / 1e3, "Kb")) +
    theme_bw() +
    coord_cartesian(
      xlim =
        c(
          0 - max(df$tlen, na.rm = TRUE) / 10,
          max(df$tlen, na.rm = TRUE) +
            max(df$tlen, na.rm = TRUE) / 10
        ),
      ylim = c(
        min(df$axis_scale_pos) - 0.25,
        max(df$axis_scale_pos) + 0.25
      )
    ) +
    scale_fill_manual(values = c(
      "Reference" = "grey",
      "Covered" = viridis_pal()(4)[3]
    ),
    na.value = "grey") +
    scale_linetype_manual(values = c(
       "Primary" = "solid",
       "Secondary" = "dashed"),
       ) +
    labs(
      x = "Position on cluster reference assembly",
      y = "Contig ID",
      fill = "Mapped\nStrand"
    ) +
    guides(linetype = "none") +
    new_scale_fill() +
    geom_rect(aes(
      xmin = qstart, xmax = qend,
      ymin = max(df$axis_scale_pos) - 0.2,
      ymax = max(df$axis_scale_pos) + 0.2,
      y = max(df$axis_scale_pos),
      fill = type,
      linetype = "solid"),
      alpha = 0.6,
      data = genes_df,
    ) +
    geom_text_repel(
      aes(y=max(df$axis_scale_pos),
          x = qstart + (qend-qstart)/2,
          xmin = qstart, xmax = qend,
          ymin = max(df$axis_scale_pos) - 0.2,
          ymax = max(df$axis_scale_pos) + 0.2,
          label = gene,
          fill = type,
          linetype = "solid",
          label_padding = 0.35,
          color = type),
      direction = "y",
      data = genes_df) +
    scale_fill_manual(
      values = c("is" = "blue",
      "amr" = "red")) +
    scale_color_manual(
      values = c("is" = "blue",
                 "amr" = "red")) +
    labs(fill = "Gene\nType",
         color = "Gene\nType")


ggplot() +
  geom_rect(aes(
    xmin = qstart, xmax = qend,
    ymin = 1 - 0.3,
    ymax = 1 + 0.3,
    fill = type,
    linetype = NA),
    data = genes_df
  )


