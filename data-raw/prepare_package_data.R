# Load raw data and save it as .rda files for blantyreESBL package

# the load_phd_data scripts are from the Thesis repo https://github.com/joelewis101/thesis

source("/Users/joelewis/Documents/PhD/Thesis/bookdown/final_cleaning_scripts/load_PhD_data.R")

source("/Users/joelewis/Documents/PhD/Thesis/bookdown/final_cleaning_scripts/load_and_clean_lims.R")

# baseline characteristics -------------------------------------------

enroll %>%
  filter(arm != 4) %>%
  mutate(
    art_time =  as.numeric(
      (data_date - hivartstart) / (365.25 /  12)
      ),
    hivart = if_else(hivart != "5A", "ART regimen: other", hivart),
    animalskept = if_else(animalskept == "Other", "Elsewhere", animalskept),
    arm = as.character(arm),
    enroll_date = data_date,
  ) %>% select(pid,
    arm,
    enroll_date,
    calc_age,
    ptsex,
    hivstatus,
    recieved_prehosp_ab,
    pmhxrechospital,
    tbstatus,
    tbongoing,
    hivonart,
    hivart,
    art_time,
    hivcpt,
    tobaccoyn,
    alcoholyn,
    highestedu,
    job,
    housholdadultsno,
    householdchildno,
    toilet ,
    watersource,
    watertreated,
    electricityyn,
    fuel,
    keepanim,
    keep.poultry,
    keep.goats,
    keep.dogs,
    keep.cattle,
    keep.sheep,
    keep.mules
  ) -> btESBL_participants


use_data(btESBL_participants, overwrite = TRUE)

### longitudinal exposure (covariate) data ---------------------------------

read_csv(
  "~/Documents/PhD/Thesis/bookdown/data/longit_covariate_data.csv") ->
  btESBL_exposures

# add in cotrim

btESBL_exposures %>%
  left_join(select(enroll, pid, hivcpt) %>%
              unique()) %>%
  mutate(cotri = if_else(hivcpt == "Yes" & !is.na(hivcpt)
                           , 1, cotri)) %>%
  select(-hivcpt) -> btESBL_exposures

use_data(btESBL_exposures, overwrite = TRUE)

# ESBL carriage data


lims_dates %>%
  rename("t" = "assess_type") ->
  btESBL_stoolESBL


use_data(btESBL_stoolESBL, overwrite = TRUE)

lims_orgs %>%
  select(lab_id,organism) %>%
  unique() %>%
  filter(!is.na(organism)) -> btESBL_stoolorgs

use_data(btESBL_stoolorgs, overwrite = TRUE)

# stan_data

read_csv("/Users/joelewis/Documents/PhD/Thesis/bookdown/data/stan_df.csv") ->
  btESBL_modeldata

btESBL_modeldata %>%
  select(-X1) -> btESBL_modeldata

use_data(btESBL_modeldata)

# stan models ----------------------------------------------

readRDS("/Users/joelewis/Documents/PhD/Thesis/bookdown/chapter_9/stan_models/model_2/stan_data_m2.rds") ->
  btESBL_model2data

use_data(btESBL_model2data, overwrite = TRUE)
"
readRDS("/Users/joelewis/Documents/PhD/Thesis/bookdown/chapter_9/stan_models/model_1/stan_data_m1.rds") ->
  btESBL_model1data

use_data(btESBL_model1data, overwrite = TRUE)

readRDS("/Users/joelewis/Documents/PhD/Thesis/bookdown/chapter_9/stan_models/model_2/stanfit_m2.rds") ->
  btESBL_model2posterior

use_data(btESBL_model2posterior, overwrite = TRUE)

readRDS("/Users/joelewis/Documents/PhD/Thesis/bookdown/chapter_9/stan_models/model_1/stanfit_m1.rds") ->
  btESBL_model1posterior

use_data(btESBL_model1posterior, overwrite = TRUE)

# simulations from posterior

btESBL_model2simulations <- read_csv("~/Documents/PhD/Thesis/bookdown/chapter_9/simulations2.csv")

use_data(btESBL_model2simulations, overwrite = TRUE)

