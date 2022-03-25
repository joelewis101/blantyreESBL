library(tidyverse)
library(blantyreESBL)

sampls <- read.csv(
  paste0(
    "~/Documents/PhD/Thesis/bookdown/",
  "chapter_5/ESBL_carriage_recruitment.csv")
  , stringsAsFactors = F)

left_join(
  sampls,
  btESBL_stoolorgs %>%
    left_join(
      select(
        btESBL_stoolESBL,
        lab_id,
        arm,
        visit
      ),
      by = "lab_id"
    ) %>%
    filter(
      ESBL == "Positive",
      # grepl("Escherichia|iella pneum", organism)
    ) %>%
    mutate(organism =
             if_else(!grepl("Escherichia|iella pneum", organism),
                     "Other",
                     organism)) %>%
    unique() %>%
    group_by(arm, visit, organism) %>%
    tally() %>%
    transmute(
      visit =
        case_when(
          visit == 1 ~ 7,
          visit == 2 ~ 28,
          visit == 3 ~ 90,
          visit == 4 ~ 180,
          visit == 0 ~ 0
        ),
      arm = paste0("Arm ", arm),
      name = case_when(
        organism == "Escherichia coli" ~ "e_coli",
        organism == "Klebsiella pneumoniae" ~ "KpSC",
        TRUE ~ organism),
      value = n
    ) %>%
    pivot_wider(id_cols = c("visit", "arm"),
                names_from = name,
                values_from = value) %>%
    filter(!is.na(visit)),
  by = c("visit","arm")) %>%

  left_join(
btESBL_sequence_sample_metadata %>%
  transmute(
      visit =
        case_when(
          visit == 1 ~ 7,
          visit == 2 ~ 28,
          visit == 3 ~ 90,
          visit == 4 ~ 180,
          visit == 0 ~ 0
        ),
      arm = paste0("Arm ", arm),
      organism =
        case_when(
          grepl("Kleb", species) ~ "KpSC.seq",
          species == "E. coli" ~ "e_coli.seq")) %>%
  group_by(organism, visit, arm) %>%
  tally() %>%
  rename(
    name = organism,
    value = n) %>%
    pivot_wider(id_cols = c("visit", "arm"),
                names_from = name,
                values_from = value)  %>%
  filter(!is.na(visit)),
  by = c("visit","arm")) %>%
  mutate(
    e_coli  = e_coli - e_coli.seq,
    KpSC = KpSC - KpSC.seq) %>%
  pivot_longer(-c("arm","visit")) %>%
  mutate(
    fill = case_when(
      grepl("seq", name) ~ "Sequenced, passed QC",
      TRUE ~ NA_character_),
    name = gsub("\\.seq","", name)) %>%
  mutate(
    name = case_when(
      name == "eligible" ~ "Participants",
      name == "collected" ~ "Samples",
      name == "e_coli" ~ "ESBL E. coli",
      name == "KpSC" ~ "ESBL KpSC",
      name == "Other" ~ "Other ESBL"),

    name = factor(name, levels =
                         c("Participants",
                           "Samples",
                           "ESBL E. coli",
                           "ESBL KpSC",
                           "Other ESBL")
  ),
  visit = factor(paste0("Day ", visit),
                 levels = paste0("Day ",
                                 c(0,7,28,90,180))
                 )) %>%
  ggplot(
    aes(fct_rev(name),value, fill = fill)) +
  geom_col() +
  coord_flip() +
  facet_grid(visit ~ arm) +
  theme_bw() +
  theme(
    text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  labs(y = "Number of participants/samples/organisms",
       x = "",
       fill = "") +
  scale_fill_manual(
    na.value = "#999999",
    values = "#39486BFF",
    breaks = "Sequenced, passed QC") -> studyflowplot

ggsave(
       here("figures/long-modelling/SF_studyflowplot.svg"),
       studyflowplot,
           width = 4, height = 4)
