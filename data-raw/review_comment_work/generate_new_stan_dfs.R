# set up dataframes

library(blantyreESBL)
library(tidyverse)
library(cmdstanr)




thesis_dir <- "~/Documents/PhD/Thesis/bookdown/"

source(paste0(thesis_dir,
              "other_scripts/stan_helpers/arrange_stan_df_functions.R"))
source(paste0(thesis_dir,
              "other_scripts/panel_data_helpers/shuffle_a_in_b2.R"))
source(paste0(thesis_dir,
              "other_scripts/panel_data_helpers/splice_ESBL2.R"))




btESBL_stoolESBL %>%
  transmute(
    pid = pid,
    lab_id = lab_id,
    sample_type = sample_type,
    ESBL = ESBL,
    assess_type = t) %>%
  group_by(pid) %>%
  do(splice_ESBL2(btESBL_exposures,. )) -> spliced



spliced$ESBL[spliced$ESBL == "Positive"] <- 1
spliced$ESBL[spliced$ESBL == "Negative"] <- 0
spliced$ESBL[is.na(spliced$ESBL)] <- 999

spliced %>%
  rowwise() %>%
  mutate(
    nonCROab = sum(
      amoxy,
      genta,
      azithro,
      tb,
      benzy,
      # cefo,
      chlora,
      cipro,
      coamo,
      clinda,
      doxy,
      erythro,
      fluclox,
      metro,
      strepto,
      cotri
    )
  ) %>% ungroup() %>%
  mutate(nonCROab = case_when(nonCROab > 0 ~ 1,
                         nonCROab == 0 ~ 0)) -> spliced



spliced  %>%
  group_by(pid) %>%
  do(generate_stan_df(., "ESBL", c("hosp",  "cefo", "nonCROab"))) -> stan.df


# zero times, code missings appropriately
stan.df[c("tstart",
          "tstop",
          grep("_time", names(stan.df), value = TRUE))] <-
  stan.df[c("tstart",
            "tstop",
            grep("_time", names(stan.df), value = TRUE))] - stan.df$tstart

stan.df %>%
  mutate(across(matches("hosp|cefo|nonCROab") &
                  matches("start_time|end_time"),
                ~ if_else(is.na(.x), -999, .x))) %>%
  mutate(across(matches("hosp|cefo|nonCROab") &
                  matches("stop_time"),
                ~ if_else(is.na(.x), 999, .x))) %>%
  select(!matches("exposure")) %>%
  ungroup() ->
  stan.df

# get to a list

N <- nrow(stan.df)
t <- stan.df$tstop

cov_mat <-
  as.matrix(dplyr::select(stan.df,
                          cefo_start_time,
                          cefo_end_time,
                          prev_cefo_stop_time,
                          nonCROab_start_time,
                          nonCROab_end_time,
                          prev_nonCROab_stop_time,
                          hosp_start_time,
                          hosp_end_time,
                          prev_hosp_stop_time))




start_state = as.matrix(
  stan.df %>%
    mutate(p0 = case_when(
      ESBL_start == 0 ~ 1,
      TRUE ~ 0),
      p1 = case_when(
        ESBL_start == 1 ~ 1,
        TRUE ~ 0)) %>%
    dplyr::select(matches("ESBL|p0|p1")) %>%
    dplyr::select(p0,p1))

end_state = as.numeric(stan.df$ESBL_stop)

n_covs <- c(0,1,2)

covs_type <- c(3,3,2)

stan_data <-
  list(
    N = N,
    t = t,
    n_covs = n_covs,
    covs_type = covs_type,
    start_state = start_state,
    end_state = end_state,
    cov_mat = covs_mat
  )

write_stan_json(
  stan_data,
  "~/.cmdstanr/cmdstan-2.28.2/models/DASSIM_CROvsall/CROvsall_data.json"
)

mod <- cmdstan_model("~/.cmdstanr/cmdstan-2.28.2/models/DASSIM_CROvsall/ESBLmod_finalV1.0_rk45.stan")

mod$sample(
  data = stan_data,
  chains = 1,
  parallel_chains = 1,
  iter_warmup = 1,
  iter_sampling = 1) -> fitz

mod$sample(
  data = btESBL_stanmodeldata,
  chains = 1,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 500) -> fitz2

file = system.file("extdata",
                   "ESBLmod_finalV1.0_rk45.stan",
                   package = "blantyreESBL")

rstan::stan(
  file = "~/.cmdstanr/cmdstan-2.28.2/models/DASSIM_CROvsall/ESBLmod_finalV1.0_rk45.stan",
  data = btESBL_stanmodeldata,
  warmup = 1,
  iter = 2,
  chains = 1) -> fitz


