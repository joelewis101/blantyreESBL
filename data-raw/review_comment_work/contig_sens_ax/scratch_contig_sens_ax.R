
# look at effect of changing cd hit settings on contig clusters

library(readr)
library(here)
library(purrr)
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)

# load data
source(
  here(paste0(
      "data-raw/review_comment_work/",
      "contig_sens_ax/scratch_load_contig_sens_ax_data.R" )
      )
  )

# load original

do.call(
  bind_rows,
  map(
    list.files(
      here("data-raw/contig_clusters/"),
      pattern = "tsv",
      full.names = TRUE
    ),
    ~ read_tsv(.x,
      col_types = "ccccccc",
      id = "file"
    )
  )
) -> df_orig

df_orig %>%
  mutate(gene = str_extract(
    file,
    "(?<=contigs-).+(?=\\.clust)"
  )) -> df_orig

# plotting functions ----------------------

add_next_edge_layer <- function(edges, dfb, dfc) {
  left_join(
    dfb %>%
      filter(clstr %in% edges$to) %>%
      select(id, clstr) %>%
      group_by(clstr) %>%
      mutate(clst.size = n()),
    dfc %>%
      select(id, clstr) %>%
      group_by(clstr) %>%
      mutate(clst.size = n()),
    by = "id"
  ) %>%
    group_by(clstr.x, clst.size.x, clstr.y, clst.size.y) %>%
    tally() %>%
    transmute(
      from_clst_size = clst.size.x,
      from = clstr.x,
      to_clst_size = clst.size.y,
      to = clstr.y,
      n = n
    ) ->
  edges_new

  return(edges_new)
}

add_new_nodes <- function(edges, dfc) {
  dfc %>%
    filter(clstr %in% edges$to) %>%
    group_by(clstr, len_diff_cutoff, ident_cutoff) %>%
    summarise(size = n()) -> nodes_new
  return(nodes_new)
}


return_cluster_plot_df <-
  function(df_new, df_orig, gene_to_plot, cluster_cutoff,
           len_cutoff_a = 0, ident_cutoff_a = 0.95,
           len_cutoff_b = 0.2, ident_cutoff_b = 0.95,
           len_cutoff_c = 0.4, ident_cutoff_c = 0.95,
           len_cutoff_d = 0.6, ident_cutoff_d = 0.95,
           len_cutoff_e = 0.8, ident_cutoff_e = 0.95) {

    # terrible code
    # df_new is the sensitivity analysis df
    # df_orig is the original data

    df_ctxm_15 <-
      df %>%
      filter(gene == gene_to_plot) %>%
      transmute(
        id = gsub("^\\.", "", contig),
        clstr = as.character(cluster),
        gene = gsub("_", "", gene),
        len_diff_cutoff = len_diff_cutoff,
        ident_cutoff = ident_cutoff
      ) %>%
      mutate(clstr = paste0(
        clstr,
        "_",
        gsub("_", "", gene),
        "_l",
        len_diff_cutoff,
        "_i",
        ident_cutoff
      ))

    dfa <-
      df_orig %>%
      filter(gene == gene_to_plot) %>%
      select(id, clstr, gene) %>%
      mutate(
        len_diff_cutoff = len_cutoff_a,
        ident_cutoff = ident_cutoff_a,
        id = gsub("^\\.", "", id)
      ) %>%
      mutate(clstr = paste0(
        clstr,
        "_",
        gsub("_", "", gene),
        "_l",
        len_diff_cutoff,
        "_i",
        ident_cutoff
      ))

    dfb <- df_ctxm_15 %>%
      filter(
        len_diff_cutoff == len_cutoff_b,
        ident_cutoff == ident_cutoff_b
      )
    left_join(
      dfa %>%
        select(id, clstr) %>%
        group_by(clstr) %>%
        mutate(clst.size = n()) %>%
        filter(clst.size > cluster_cutoff),
      dfb %>%
        select(id, clstr) %>%
        group_by(clstr) %>%
        mutate(clst.size = n()),
      by = "id"
    ) %>%
      group_by(clstr.x, clst.size.x, clstr.y, clst.size.y) %>%
      tally() %>%
      transmute(
        from_clst_size = clst.size.x,
        from = clstr.x,
        to_clst_size = clst.size.y,
        to = clstr.y,
        n = n
      ) ->
    edges
    bind_rows(
      dfa %>%
        filter(clstr %in% edges$from),
      dfb %>%
        filter(clstr %in% edges$to)
    ) %>%
      group_by(
        clstr, len_diff_cutoff,
        ident_cutoff
      ) %>%
      summarise(size = n()) -> nodes
if (!is.na(len_cutoff_c) & !is.na(ident_cutoff_c)) {
    dfc <- df_ctxm_15 %>%
      filter(
        len_diff_cutoff == len_cutoff_c,
        ident_cutoff == ident_cutoff_c
      )
    add_next_edge_layer(edges, dfb, dfc) -> newedges
    add_new_nodes(newedges, dfc) -> newnodes
    nodes <- bind_rows(nodes, newnodes)
    edges <- bind_rows(edges, newedges)
}

if (!is.na(len_cutoff_d) & !is.na(ident_cutoff_d)) {
    dfd <- df_ctxm_15 %>%
      filter(
        len_diff_cutoff == len_cutoff_d,
        ident_cutoff == ident_cutoff_d
      )
    add_next_edge_layer(edges, dfc, dfd) -> newedges
    add_new_nodes(newedges, dfd) -> newnodes
    nodes <- bind_rows(nodes, newnodes)
    edges <- bind_rows(edges, newedges)
}

if (!is.na(len_cutoff_e) & !is.na(ident_cutoff_e)) {
    dfe <- df_ctxm_15 %>%
      filter(
        len_diff_cutoff == len_cutoff_e,
        ident_cutoff == ident_cutoff_e
      )
    add_next_edge_layer(edges, dfd, dfe) -> newedges
    add_new_nodes(newedges, dfe) -> newnodes
    nodes <- bind_rows(nodes, newnodes)
    edges <- bind_rows(edges, newedges)
}
    tbl_graph(
      nodes = nodes,
      edges = edges,
      directed = TRUE
    ) -> ctxm15_g
    return(ctxm15_g)
  }

# start plotting bruv ------------------------------

# cluster cumulative membership ----------------

df %>%
  filter( ident_cutoff == 0.95) %>%
  transmute(
    id = gsub("^\\.","", contig),
    clstr = as.character(cluster),
    gene = gsub("_","",gene),
    len_diff_cutoff = len_diff_cutoff,
    ident_cutoff = ident_cutoff) %>%
  mutate(clstr = paste0(
    clstr,
    "_",
    gsub("_","",gene),
    "_l",
    len_diff_cutoff,
    "_i",
    ident_cutoff
  ) ) %>%
  bind_rows(
     df_orig %>%
              select(id, clstr, gene) %>%
              mutate(len_diff_cutoff = 0,
                     ident_cutoff = 0.95,
                     id = gsub("^\\.","", id),
                     gene = gsub("_","", gene)) %>%
              mutate(clstr = paste0(
                clstr,
                "_",
                gsub("_","",gene),
                "_l",
                len_diff_cutoff,
                "_i",
                ident_cutoff
              ))) %>%
  group_by(len_diff_cutoff, gene) %>%
  mutate(
    n_gene = n()) %>%
  group_by(gene) %>%
    mutate(
      gene = paste0(gene, " (n = ", max(n_gene),")")) %>%
  ungroup() %>%
  mutate(
    gene = fct_infreq(gene)) %>%
  group_by(len_diff_cutoff, clstr, gene) %>%
  summarise(n = n()) %>%
  arrange(gene,len_diff_cutoff, -n) %>%
  group_by(gene,len_diff_cutoff) %>%
  mutate(n_clust = 1:n()) %>%
  mutate(
    prop = cumsum(n)/sum(n)) %>%
  ggplot(aes(n_clust,
             prop,
             group = len_diff_cutoff,
             color = len_diff_cutoff)) +
  geom_step() +
  facet_wrap(~ gene, scales = "free_x")

df %>%
  filter(len_diff_cutoff == 0) %>%
  transmute(
    id = gsub("^\\.","", contig),
    clstr = as.character(cluster),
    gene = gsub("_","",gene),
    len_diff_cutoff = len_diff_cutoff,
    ident_cutoff = ident_cutoff) %>%
  mutate(clstr = paste0(
    clstr,
    "_",
    gsub("_","",gene),
    "_l",
    len_diff_cutoff,
    "_i",
    ident_cutoff
  ) ) %>%
  bind_rows(
     df_orig %>%
              select(id, clstr, gene) %>%
              mutate(len_diff_cutoff = 0,
                     ident_cutoff = 0.95,
                     id = gsub("^\\.","", id),
                     gene = gsub("_","", gene)) %>%
              mutate(clstr = paste0(
                clstr,
                "_",
                gsub("_","",gene),
                "_l",
                len_diff_cutoff,
                "_i",
                ident_cutoff
              ))) %>%
  group_by(ident_cutoff, gene) %>%
  mutate(
    n_gene = n()) %>%
  group_by(gene) %>%
    mutate(
      gene = paste0(gene, " (n = ", max(n_gene),")")) %>%
  ungroup() %>%
  mutate(
    gene = fct_infreq(gene)) %>%
  group_by(ident_cutoff, clstr, gene) %>%
  summarise(n = n()) %>%
  arrange(gene,ident_cutoff, -n) %>%
  group_by(gene,ident_cutoff) %>%
  mutate(n_clust = 1:n()) %>%
  mutate(
    prop = cumsum(n)/sum(n)) %>%
  ggplot(aes(n_clust,
             prop,
             group = ident_cutoff,
             color = ident_cutoff)) +
  geom_step() +
  facet_wrap(~ gene, scales = "free_x")

# plot cluster sizes

return_cluster_plot_df(
  df,
  df_orig,
  "CTX_M_27",
  5
) %>%
  ggraph(layout = "sugiyama") +
  geom_edge_link(aes(edge_width = n), alpha = 0.7) +
  geom_node_point(aes(size = size, color = as.factor(len_diff_cutoff)))



  return_cluster_plot_df(
  df,
  df_orig,
  "CTX_M_15",
  5,
  len_cutoff_a = 0, ident_cutoff_a = 0.95,
  len_cutoff_b = 0, ident_cutoff_b = 0.96,
  len_cutoff_c = 0, ident_cutoff_c = 0.97,
  len_cutoff_d = 0, ident_cutoff_d = 0.99,
  len_cutoff_e = 0, ident_cutoff_e = 1
)  %>%
  ggraph(layout = "sugiyama") +
  geom_edge_link(aes(edge_width = n), alpha = 0.7) +
  geom_node_point(aes(size = size, color = as.factor(ident_cutoff)))

return_cluster_plot_df(
  df,
  df_orig,
  "CTX_M_27",
  5
) %>%
  ggraph(layout = "sugiyama") +
  geom_edge_link(aes(edge_width = n), alpha = 0.7) +
  geom_node_point(aes(size = size, color = as.factor(len_diff_cutoff)))

return_cluster_plot_df(
  df,
  df_orig,
  "CTX_M_27",
  5,
  len_cutoff_a = 0, ident_cutoff_a = 0.95,
  len_cutoff_b = 0, ident_cutoff_b = 0.96,
  len_cutoff_c = 0, ident_cutoff_c = 0.97,
  len_cutoff_d = 0, ident_cutoff_d = 0.99,
  len_cutoff_e = 0, ident_cutoff_e = 1
)  %>%
  ggraph(layout = "sugiyama") +
  geom_edge_link(aes(edge_width = n), alpha = 0.7) +
  geom_node_point(aes(size = size, color = as.factor(ident_cutoff)))
