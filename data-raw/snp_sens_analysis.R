# sensitivity analysis of snp cutoffs

btESBL_snpdists_esco %>%
  pivot_longer(-sample) %>%
  rename("sample.x" = "sample",
         "sample.y" = "name",
         "snpdist_esco" = "value") %>%
  left_join(
    select(btESBL_sequence_sample_metadata, lane, supplier_name) %>%
      rename("lab_id.x" = "supplier_name"),
    by = c("sample.x" = "lane")) %>%
  left_join(
    select(btESBL_sequence_sample_metadata, lane, supplier_name) %>%
      rename("lab_id.y" = "supplier_name"),
    by = c("sample.y" = "lane")) %>%
  group_by(lab_id.x, lab_id.y) %>%
  slice(n=1) -> snpdist.e.long



btESBL_snpdists_kleb  %>%
  pivot_longer(-sample) %>%
  rename("sample.x" = "sample",
         "sample.y" = "name",
         "snpdist_kleb" = "value") %>%
  left_join(
    select(btESBL_sequence_sample_metadata, lane, supplier_name) %>%
      rename("lab_id.x" = "supplier_name"),
    by = c("sample.x" = "lane")) %>%
  left_join(
    select(btESBL_sequence_sample_metadata, lane, supplier_name) %>%
      rename("lab_id.y" = "supplier_name"),
    by = c("sample.y" = "lane")) %>%
  group_by(lab_id.x, lab_id.y) %>%
  slice(n=1)  -> snpdist.k.long

# make list-column df -----------------------------------------------------


btESBL_stoolESBL %>%
  left_join(select(btESBL_stoolorgs, lab_id, organism),
            by = "lab_id") %>%
  nest(orgs = organism) %>%
  left_join(
    btESBL_popPUNK %>%
      left_join(
        select(btESBL_sequence_sample_metadata , lane, supplier_name),
        by = c("Taxon" = "lane")
      ) %>%
      select(supplier_name, Cluster),
    by = c("lab_id" = "supplier_name")
  ) %>%
  nest(pp_clust = Cluster) %>%
  left_join(
    btESBL_contigclusters %>%
      left_join(
        select(btESBL_sequence_sample_metadata , lane, supplier_name)
      ) %>%
      select(supplier_name, clstr_name),
    by = c("lab_id" = "supplier_name")
  ) %>%
  nest(contig_clust = clstr_name) %>%
  left_join(
    select(btESBL_sequence_sample_metadata, lane, supplier_name),
    by = c("lab_id" = "supplier_name")
  ) %>%
  nest(lanes = lane) -> samples

samples %>%
  full_join(samples, by = character())  %>%
  filter(lab_id.x != lab_id.y) %>%
  mutate(delta_t =
           interval(data_date.x, data_date.y) / days(1)) %>%
  filter(delta_t >= 0) -> samples
# compare presence absence of clusters

samples %>%
  mutate(esbl.x.and.y = ESBL.x == "Positive" &
           ESBL.y == "Positive") %>%
  #compare ESCO poppunk clusters between x and y
  # and add variables for e coli cluster and presence/absence
  # e coli to later remove isolates that weren't sequenced
  mutate(
    same.esco.poppunk.xandy =
      map2(pp_clust.x,
           pp_clust.y,
           ~ .x$Cluster[grepl("E", .x$Cluster)] %in%
             .y$Cluster[grepl("E", .y$Cluster)]) %>%
      map_lgl(any),
    # flags for existence of poppunk clusters
    esco.poppunk.cluster.exists.x =
      map(pp_clust.x, ~ grepl("E", .x$Cluster)) %>%
      map_lgl(any),
    esco.poppunk.cluster.exists.y =
      map(pp_clust.y, ~ grepl("E", .x$Cluster)) %>%
      map_lgl(any),
    # flag for existence of esco
    esco.exists.x =
      map_lgl(orgs.x, ~ any(grepl("coli", .x$organism))),
    esco.exists.y =
      map_lgl(orgs.y,  ~ any(grepl("coli", .x$organism))),
    same.esco.xandy = esco.exists.x & esco.exists.y,
    # contig clusters
    same.contig.cluster =
      map2(contig_clust.x,
           contig_clust.y,
           ~ .x$clstr_name %in%
             .y$clstr_name) %>%
      map_lgl(any)
  ) %>%
  # same but for klebs
  mutate(
    same.kleb.poppunk.xandy =
      map2(pp_clust.x,
           pp_clust.y,
           ~ .x$Cluster[grepl("K", .x$Cluster)] %in%
             .y$Cluster[grepl("K", .y$Cluster)]) %>%
      map_lgl(any),
    # flags for existence of poppunk clusters
    kleb.poppunk.cluster.exists.x =
      map(pp_clust.x, ~ grepl("K", .x$Cluster)) %>%
      map_lgl(any),
    kleb.poppunk.cluster.exists.y =
      map(pp_clust.y, ~ grepl("K", .x$Cluster)) %>%
      map_lgl(any),
    # flag for existence of kleb
    kleb.exists.x =
      map_lgl(orgs.x, ~ any(grepl("Klebsiella pneumoniae", .x$organism))),
    kleb.exists.y =
      map_lgl(orgs.y,  ~ any(grepl("Klebsiella pneumoniae", .x$organism))),
    same.kleb.xandy = kleb.exists.x & kleb.exists.y) ->
  samples

# merge in snpdist clusters

samples %>%
  left_join(
    select(snpdist.e.long, lab_id.x, lab_id.y, snpdist_esco),
    by = c("lab_id.x", "lab_id.y")
  ) %>%
  left_join(
    select(snpdist.k.long, lab_id.x, lab_id.y, snpdist_kleb),
    by = c("lab_id.x", "lab_id.y")
  ) -> samples

samples%>%
  filter(pid.x == pid.y)-> pairwise_within


pairwise_within %>%
  filter(esco.exists.x) %>%
  select(data_date.x, data_date.y,
         delta_t,
         esco.exists.x,
         esco.poppunk.cluster.exists.x,
         esco.exists.y,
         esco.poppunk.cluster.exists.y,
         delta_t,
         snpdist_esco) %>%
  pivot_longer(-c(data_date.x, data_date.y,
                  esco.exists.x, esco.exists.y,
                  esco.poppunk.cluster.exists.x,
                  esco.poppunk.cluster.exists.y,
                  delta_t)) %>%
  #filter out those with an esco but no poppunk cluster - they've not been
  # sequenced
  filter(
    !(esco.exists.x &
      !esco.poppunk.cluster.exists.x),
    !(esco.exists.y &
      !esco.poppunk.cluster.exists.y)
  ) -> ecoli.long

for (i in 0:20) {
  newcolvar <- paste0("snpdist_",i)
  ecoli.long %>%
    mutate({{newcolvar}} := if_else(
      value <= i & !is.na(value), TRUE, FALSE)
    ) -> ecoli.long
}

ecoli.long %>%
  select(delta_t, contains("snpdist")) %>%
  pivot_longer(-delta_t) %>%
  mutate(snpdist = as.numeric(gsub("snpdist_","", name))) %>%
  ggplot(aes(delta_t, as.numeric(value), color = snpdist, group = name)) +
    geom_smooth(se = FALSE) +
  scale_color_viridis_c() +
  coord_cartesian(xlim = c(0,150), ylim = c(0,0.13)) +
  theme_bw() +
  scale_y_continuous(breaks = c(0,0.02,0.04,0.06,0.08,0.1,0.12)) +
  labs(color = "SNP\ndistance",
       x = "Time/days",
       y = "Proportion") -> a



## klebs


pairwise_within %>%
  filter(kleb.exists.x) %>%
  select(data_date.x, data_date.y,
         delta_t,
         kleb.exists.x,
         kleb.poppunk.cluster.exists.x,
         kleb.exists.y,
         kleb.poppunk.cluster.exists.y,
         delta_t,
         snpdist_kleb) %>%
  pivot_longer(-c(data_date.x, data_date.y,
                  kleb.exists.x, kleb.exists.y,
                  kleb.poppunk.cluster.exists.x,
                  kleb.poppunk.cluster.exists.y,
                  delta_t)) %>%
  #filter out those with an kleb but no poppunk cluster - they've not been
  # sequenced
  filter(
    !(kleb.exists.x &
        !kleb.poppunk.cluster.exists.x),
    !(kleb.exists.y &
        !kleb.poppunk.cluster.exists.y)
  ) -> kleb.long

for (i in 0:20) {
  newcolvar <- paste0("snpdist_",i)
  kleb.long %>%
    mutate({{newcolvar}} := if_else(
      value <= i & !is.na(value), TRUE, FALSE)
    ) -> kleb.long
}

kleb.long %>%
  select(delta_t, contains("snpdist")) %>%
  pivot_longer(-delta_t) %>%
  mutate(snpdist = as.numeric(gsub("snpdist_","", name))) %>%
  ggplot(aes(delta_t, as.numeric(value), color = snpdist, group = name)) +
  geom_smooth(se = FALSE) +
  scale_color_viridis_c() +
  coord_cartesian(xlim = c(0,150),ylim = c(0,0.13)) +
  theme_bw() +
  scale_y_continuous(breaks = c(0,0.02,0.04,0.06,0.08,0.1,0.12)) +
  labs(color = "SNP\ndistance",
       x = "Time/days",
       y = "Proportion") -> b

a + b + plot_annotation(tag_levels = "A") + plot_layout(guides = "collect")

## snp network

make_cluster_pairwise_comparison_df <- function(pairwise_snp_df,
                                                metadata_df,
                                                cut_tree_vect,
                                                n) {
  cluster_samples <- names(cut_tree_vect[cut_tree_vect == n])
  pairwise_snp_df <- as.data.frame(pairwise_snp_df)
  rownames(pairwise_snp_df) <- pairwise_snp_df$sample
  pairwise_snp_df[cluster_samples, c("sample", cluster_samples)] %>%
    pivot_longer(-sample) %>%
    filter(sample != name) -> d

  d[!duplicated(apply(d[1:2], 1, sort), MARGIN = 2), ] -> d
  names(d) <- c("sample.x", "sample.y", "snpdist")
  d$sample.x <- gsub("#", "_", d$sample.x)
  d$sample.y <- gsub("#", "_", d$sample.y)


  left_join(
    d,
    select(
      metadata_df,
      lane,
      pid,
      arm,
      visit,
      enroll_date,
      hospoutcomedate,
      data_date
    ),
    by = c("sample.x" = "lane")
  ) %>%
    left_join(
      select(
        metadata_df,
        lane,
        pid,
        arm,
        visit,
        enroll_date,
        hospoutcomedate,
        data_date
      ),
      by = c("sample.y" = "lane")
    ) -> d
  d$pair <- paste0(d$sample.x, "-", d$sample.y)
  d$cluster_number <- n
  d$cluster_name <- n
  d %>%
    mutate(type = case_when(pid.x == pid.y ~ "within",
                            TRUE ~ "between")) ->   d
  return(d)

}


make_all_clusters_pairwise_comparison_df  <-
  function(pairwise_snp_df,
           metadata_df,
           cut_tree_vect) {
    out <- list()
    for (i in 1:max(cut_tree_vect)) {
      #    print(i)
      #  print(i)
      make_cluster_pairwise_comparison_df(pairwise_snp_df,
                                          metadata_df,
                                          cut_tree_vect,
                                          i) -> out[[i]]

    }
    return(do.call(rbind, out))


  }

#### start plots -------------

outlist.e.edges <- list()
outlist.e.vertices <- list()
outlist.e.plots <- list()

outlist.k.edges <- list()
outlist.k.vertices <- list()
outlist.k.plots <- list()
listindex <- 0
for (i in c(0,3,7,10)) {
  listindex = listindex + 1
  hclust(as.dist(btESBL_snpdists_esco[-1])) -> hclust_snpdists.e
  cutree(hclust_snpdists.e, h = i) -> cut_tree_vect.e


  hclust(as.dist(btESBL_snpdists_kleb[-1])) -> hclust_snpdists.k
  cutree(hclust_snpdists.k, h = i) -> cut_tree_vect.k

  print(i)
  as.data.frame(cut_tree_vect.e) %>%
    group_by(cut_tree_vect.e) %>%
    mutate(n = n()) %>%
    filter(n > 1) %>%
    ungroup() %>%
    summarise(med = median(n),
              lq = quantile(n, 0.25),
              uq = quantile(n, 0.75))

  as.data.frame(cut_tree_vect.k) %>%
    group_by(cut_tree_vect.k) %>%
    mutate(n = n()) %>%
    filter(n > 1) %>%
    ungroup() %>%
    summarise(med = median(n),
              lq = quantile(n, 0.25),
              uq = quantile(n, 0.75))

  make_all_clusters_pairwise_comparison_df(btESBL_snpdists_esco,
                                           btESBL_sequence_sample_metadata,
                                           cut_tree_vect.e) -> pairwise_snpclust.e


  bind_rows(
    btESBL_sequence_sample_metadata %>%
      semi_join(btESBL_snpdists_esco,
                by = c("lane" = "sample")) %>%
      select(lane, pid) %>%
      full_join(
        select(btESBL_sequence_sample_metadata, lane, pid) %>%
          semi_join(btESBL_snpdists_esco ,
                    by = c("lane" = "sample")),
        by = character()
      )  %>%
      filter(lane.x != lane.y,
             pid.x == pid.y) %>%
      select(lane.x, lane.y) %>%
      mutate(type = "Same Patient") %>%
      rename(from = lane.x,
             to = lane.y),
    pairwise_snpclust.e %>%
      select(sample.x, sample.y) %>%
      rename(from = sample.x,
             to = sample.y) %>%
      mutate(type = "SNP cluster")

  ) -> edges.e




  btESBL_sequence_sample_metadata %>%
    semi_join(btESBL_snpdists_esco,
              by = c("lane" = "sample")) %>%
    group_by(lane) %>%
    slice(n = 1) %>%
    filter(!is.na(hosp_assoc)) %>%
    mutate(
      hosp_assoc = case_when(
        hosp_assoc == "community" ~ "Community",
        hosp_assoc == "in_hospital" ~ "In hospital",
        hosp_assoc == "recent_dc" ~ "Post discharge"
      ),
    ) %>%
    relocate(lane, .before = everything()) ->
    vertices.e

  edges.e %>%
    mutate(sortvar = map2_chr(from, to, ~ paste(sort(c(
      .x, .y
    )), collapse = ""))) %>%
    group_by(sortvar, type) %>%
    slice(n = 1) %>%
    ungroup() %>%
    select(-sortvar) %>%
    group_by(from, to) %>%
    mutate(weight = 0.1) %>%
    # remove those with missing metadata
    anti_join(btESBL_sequence_sample_metadata %>%
                filter(is.na(hosp_assoc)),
              by = c("from" = "lane")) %>%
    anti_join(btESBL_sequence_sample_metadata %>%
                filter(is.na(hosp_assoc)),
              by = c("to" = "lane")) -> edges.e


  gr.e <-
    graph_from_data_frame(edges.e, directed = FALSE, vertices = vertices.e)

  ggraph(gr.e, layout = "fr") + geom_edge_fan(aes(color = type), width =
                                                1) +
    geom_node_point(aes(fill = hosp_assoc), shape = 21, size = 3) +
    #  facet_edges(~ type) +
    theme_void() + scale_fill_manual(values = c("white", "black", "grey")) +
    theme(legend.title = element_blank()) ->
    snp_network_plot.esco

  outlist.e.edges[[listindex]] <- edges.e
  outlist.e.vertices[[listindex]] <- vertices.e
  outlist.e.plots[[listindex]] <- snp_network_plot.esco



  # kleb

  make_all_clusters_pairwise_comparison_df(btESBL_snpdists_kleb,
                                           btESBL_sequence_sample_metadata,
                                           cut_tree_vect.k) -> pairwise_snpclust.k


  bind_rows(
    btESBL_sequence_sample_metadata %>%
      semi_join(btESBL_snpdists_kleb,
                by = c("lane" = "sample")) %>%
      select(lane, pid) %>%
      full_join(
        select(btESBL_sequence_sample_metadata, lane, pid) %>%
          semi_join(btESBL_snpdists_kleb ,
                    by = c("lane" = "sample")),
        by = character()
      )  %>%
      filter(lane.x != lane.y,
             pid.x == pid.y) %>%
      select(lane.x, lane.y) %>%
      mutate(type = "Same Patient") %>%
      rename(from = lane.x,
             to = lane.y),
    pairwise_snpclust.k %>%
      select(sample.x, sample.y) %>%
      rename(from = sample.x,
             to = sample.y) %>%
      mutate(type = "SNP cluster")

  ) -> edges.k




  btESBL_sequence_sample_metadata %>%
    semi_join(btESBL_snpdists_kleb,
              by = c("lane" = "sample")) %>%
    group_by(lane) %>%
    slice(n = 1) %>%
    filter(!is.na(hosp_assoc)) %>%
    mutate(
      hosp_assoc = case_when(
        hosp_assoc == "community" ~ "Community",
        hosp_assoc == "in_hospital" ~ "In hospital",
        hosp_assoc == "recent_dc" ~ "Post discharge"
      ),
    ) %>%
    relocate(lane, .before = everything()) ->
    vertices.k

  edges.k %>%
    mutate(sortvar = map2_chr(from, to, ~ paste(sort(c(
      .x, .y
    )), collapse = ""))) %>%
    group_by(sortvar, type) %>%
    slice(n = 1) %>%
    ungroup() %>%
    select(-sortvar) %>%
    group_by(from, to) %>%
    mutate(weight = n() / 10) %>%
    # remove those with missing metadata
    semi_join(vertices.k,
              by = c("from" = "lane")) %>%
    semi_join(vertices.k, by = c("to" = "lane")) -> edges.k


  gr.k <-
    graph_from_data_frame(edges.k, directed = FALSE, vertices = vertices.k)

  ggraph(gr.k, layout = "fr") + geom_edge_fan(aes(color = type), width =
                                                1) +
    geom_node_point(aes(fill = hosp_assoc), shape = 21, size = 3) +
    #  facet_edges(~ type) +
    theme_void() + scale_fill_manual(values = c("white", "black", "grey")) +
    theme(legend.title = element_blank())  ->
    snp_network_plot.kleb

  outlist.k.edges[[listindex]] <- edges.k
  outlist.k.vertices[[listindex]] <- vertices.k
  outlist.k.plots[[listindex]] <- snp_network_plot.kleb

}

outlist.e.plots[[1]] +
  outlist.e.plots[[2]] +
  outlist.e.plots[[3]] +
  outlist.e.plots[[4]] +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")

outlist.k.plots[[1]] +
  outlist.k.plots[[2]] +
  outlist.k.plots[[3]] +
  outlist.k.plots[[4]] +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")

