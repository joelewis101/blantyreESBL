# Plot paf file alignments from minimap
# JL Jan 2022


# PopGenome::readData("data-raw/review_comment_work/contig_msa/alignments/CTXM27/",
#                     format = "fasta",
#                     include.unknown = TRUE) -> msa

# This function takes
#
# 1) A tibble with alignment in paf format
# 2) A multiple sequence alignment fasta in GENOME format from
#    PopGenome
#
# And returns a ggplot with nucleotide diversity and the
# queries shown mapped to target with colours indicating coverage
# - the colour shows whether + or - mapped (as defined by PAF file)
#
# The paf tibble should have colnames as per pafR package
#     c("qname","qlen","qstart","qend","strand","tname", "tlen", "tstart",
#      "tend","nmatch","alen","mapq")
# See
# https://cran.r-project.org/
# web/packages/pafr/vignettes/Introduction_to_pafr.html
#
# plot title is just string to add at top of plot
# window_val = width of sliding window in bases
# jump_val = number of bases sliding window moves forward for next window

return_contig_cluster_plots <- function(paf_df,
                                        msa,
                                        plot_title,
                                        window_val = 1,
                                        jump_val = 1) {
  # Add start and end query locations on target
  # Need to flip if negative strand

  require(dplyr)
  require(PopGenome)
  require(ggplot2)
  require(stringr)

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

  # plot alignments
  ggplot(df, aes(
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
      "+" = viridis_pal()(4)[3],
      "-" = viridis_pal()(4)[2]
    ),
    na.value = "grey") +
    scale_linetype_manual(values = c(
      "Primary" = "solid",
      "Secondary" = "dashed"),
      na.value = "grey"
      ) +
    # theme(legend.position = "bottom") +
    labs(
      x = "Position on cluster reference assembly",
      y = "Contig ID",
      fill = "Mapped\nStrand"
    ) +
    guides(linetype = "none")  -> p

  # now nucleotide diversity
  sliding.window.transform(msa,
                           width = window_val,
                           jump = jump_val,
                           type = 2) -> shv2
  diversity.stats(shv2) -> msa2
  data.frame(msa2@nuc.diversity.within) -> nuc_div
  nuc_div$window <- rownames(nuc_div)

  # get window location for plotting
  nuc_div %>%
    mutate(
      window_left = as.numeric(str_extract(window, "^[0-9]+(?= )")),
      window_right = as.numeric(str_extract(window, "(?<=- )[0-9]+(?= )")),
      window_mid = window_left +
        (window_right - window_left) / 2
    ) %>%
    dplyr::rename(diversity = pop.1) %>%
    mutate(diversity = diversity / window_val) -> nuc_div

  # plot nuc diversity
  nuc_div %>%
    ggplot(aes(window_mid, diversity)) +
    geom_line() +
    theme_bw() +
    coord_cartesian(
      xlim =
        c(
          0 - max(df$tlen, na.rm = TRUE) / 10,
          max(df$tlen, na.rm = TRUE) +
            max(df$tlen, na.rm = TRUE) / 10
        )
    ) +
    labs(
      x = element_blank(),
      y = "Cluster\nnucleotide diversity",
      title = {{ plot_title }}
    ) +
    theme(axis.text.x = element_blank()) +
    ylim(c(0,1)) -> p2
  plot_out <- (p2 + p) + plot_layout(ncol = 1, heights = c(1, 3))
  return(plot_out)
}



