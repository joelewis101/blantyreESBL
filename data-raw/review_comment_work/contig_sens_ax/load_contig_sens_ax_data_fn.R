
# loads contig sens ax data
# load dat data ---------------------------------------

load_contig_sens_ax_data <- function() {

do.call(
  bind_rows,
  map(
    list.files(
      here("data-raw/review_comment_work/contig_sens_ax"),
      pattern = "tsv",
      full.names = TRUE
    ),
    read_tsv,
    show_col_types = FALSE
  )
) -> df

df %>%
  group_by(dir, filename) %>%
  tally() %>%
  as.data.frame()

df %>%
  transmute(
    cluster = cluster,
    id = id,
    length = length,
    contig = contig,
    lane = Lane,
    len_diff_cutoff =
      case_when(
        grepl("len_diff", dir) ~
          as.numeric(
            str_extract(dir, "(?<=len_diff_).{3}")
          ),
        TRUE ~ 0
      ),
    ident_cutoff = case_when(
      grepl("seq_ident", dir) ~
        as.numeric(
          str_extract(dir, "(?<=seq_ident_)[0-9]\\.[0-9]{1,3}")
        ),
      TRUE ~ 0.95
    ),
    gene = gsub("contigs-|\\.clust\\.clstr", "", filename)
  ) -> df

return(df)
}
