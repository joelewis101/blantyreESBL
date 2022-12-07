
library(tidyverse)
library(here)

read_tsv(here("data-raw/ecoli-genomics-paper/horesh_full_diversity/FINAL_METADATA_CLEANED.tsv")) -> df

df %>% 
  separate(Assembly_Location,into = c("Assembly_Location_1", 
                                        "Assembly_Location_2",
                                        "Assembly_Location_3"
                                        ),
           sep = ",") 


df %>%
  select(ID, Name, contains("Assembly")) %>%
  pivot_longer(-c(ID, Name)) %>%
  filter(!is.na(value)) %>%
  pull(value) %>%
  write_lines("~/tmp/horesh_assemblies.txt")

df %>%
  separate(Run_ID, into = c("run1",
                            "run2",
                            "run3",
                            "run4",
                            "run5",
                            "run6",
                            "run7",
                            "run8",
                            "run9",
                            "run10"), sep = ",") %>%
select(ID, Name, contains("run")) %>%
  pivot_longer(-c(ID, Name)) %>%
  filter(!is.na(value)) %>%
  pull(value) %>%
  write_lines("~/tmp/horesh_run_ids.txt")



df %>%
  separate(Run_ID, into = c("run"), sep = ",") %>%
select(ID, Name, contains("run"))


  pivot_longer(-c(ID, Name)) %>%
  filter(!is.na(value)) %>%
  pull(value) %>%
  write_lines("~/tmp/horesh_run_ids.txt")



df %>%
  separate(Reads_Location, into = c("reads1",
                                    "reads2"),
                                    sep = ",") %>%
filter(Make_artificial == "Yes") %>%
mutate(reads1 = paste0(gsub("_[0-9].fastq","",reads1),"*")) %>%
pull(reads1) %>%
write_lines("~/tmp/horesh_artifical_reads_filepaths.txt")

df %>%
  separate(Run_ID, into = c("run1"),sep = ",") %>%
filter(Make_artificial == "No") %>%
mutate(run = if_else(
                        grepl("trimmed",Reads_Location),
                        paste0(run1, "_trimmed"),
                        run1)) %>%
pull(run)



read_lines("~/tmp/sampleid_assembly_paths.txt") -> sampl_661k

sampl_661k[grepl("SAMN03492703", sampl_661k)]


df %>%
  separate(Reads_Location,
    into = paste0(1:4,
      "reads"),
    sep = ",") %>%
  select(ID, contains("reads")) %>%
  pivot_longer(-ID) %>%
  arrange(ID, value) %>%
  filter(!is.na(value)) %>%
  pull(value) %>% 
  write_lines(here("data-raw/ecoli-genomics-paper/horesh_full_diversity/horesh_read_filepaths.txt"))


df %>%
  separate(Assembly_Location,
    into = c(
      "Assembly_Location_1",
      "Assembly_Location_2",
      "Assembly_Location_3"
    ),
    sep = ","
  ) %>%
  select(ID, starts_with("Assembly")) %>%
  pivot_longer(-ID) %>%
  filter(!is.na(value)) %>%
  pull(value) %>%
  write_lines("~/tmp/horesh_assembly_filepaths.txt")



df %>%
  separate(Assembly_Location,
    into = c(
      "Assembly_Location_1",
      "Assembly_Location_2",
      "Assembly_Location_3"
    ),
    sep = ","
  ) %>%
  select(ID, starts_with("Assembly")) %>%
  pivot_longer(-ID) %>%
  filter(!is.na(value)) %>%
  pull(value) %>%
  write_lines("~/tmp/horesh_assembly_filepaths.txt")


# prep lanes for ariba


df %>%
  separate(Reads_Location,
    into = paste0(1:4,
      "reads"),
    sep = ",") %>%
  select(ID, contains("reads")) %>%
  pivot_longer(-ID) %>%
  group_by(ID) %>%
  arrange(ID, value) %>%
  slice(1:2) %>%
  filter(!is.na(value)) %>%
  group_by(ID) %>%
  filter(n() == 2) %>%
  pull(value) %>% 
  write_lines(here("data-raw/ecoli-genomics-paper/horesh_full_diversity/read_filepaths_for_ariba.txt"))



df %>%
  separate(Reads_Location,
    into = paste0(1:4,
      "reads"),
    sep = ",") %>%
  select(ID, contains("reads")) %>%
  pivot_longer(-ID) %>%
  group_by(ID) %>%
  arrange(ID, value) %>%
  slice(1:2) %>%
  filter(!is.na(value)) %>%
  group_by(ID) %>%
  filter(n() == 2) %>%
  filter(grepl("fastq$", value))
