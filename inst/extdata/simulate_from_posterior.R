library(rstan)
library(blantyreESBL)
library(deSolve)
library(dplyr)
#needed to save
#library(here)

## simulations of different scenarios from the posterior of the fitted models

fit_mod2 <- btESBL_model2posterior


# Make stan functions available to take advantage of the C++ code
# to quickly calculate the values of the time varying covariates
# that's the function return_time_varying_coefs_exp_flat()

expose_stan_functions(file = system.file("extdata",
                                         "ESBLmod_finalV1.0_rk45.stan",
                                         package = "blantyreESBL"),
                      fit_mod2)

# function for deSolve to solve differential eqs

ode_finalstateprob <- function(t, state, parameters) {
  coefs <-
    return_time_varying_coefs_exp_flat(parameters$cov_mat,
                                        t,
                                        parameters$covs_type,
                                        parameters$gammas)

  dp0 <-
    -(parameters$lambda * state[[1]] * exp(sum(parameters$betas * coefs)))   +
    (parameters$mu * state[[2]] * exp(sum(parameters$alphas * coefs)))
  dp1 <-
    (parameters$lambda * state[[1]] * exp(sum(parameters$betas * coefs)))  -
    (parameters$mu * state[[2]] * exp(sum(parameters$alphas * coefs)))
  return(list(c(dp0, dp1)))
}

# function to solve differential state equations for each posterior draw

iterate_over_posterior_parms <- function(p.params, t, cov_mat, states, pid) {
  purrr::pmap(p.params, ~ode(y = states, t = t, func = ode_finalstateprob,
                      parms = list(cov_mat = as.matrix(cov_mat), covs_type = c(3,2),
                                   alphas = c(..1 ,..2),
                                   betas = c(..3,..4),
                                   gammas = ..5, lambda =..6,
                                   mu =..7))) -> out
  do.call(rbind, out) -> out
  as.data.frame(out) -> out
  out$draw <- rep(1:nrow(p.params), each = length(t))
  out$pid <- pid
  return(out)
}


day_cut <- 1

# get the posterior parameters to use from the fitted model

rstan::extract(fit_mod2,
               pars = c("alphas", "betas", "gammas", "lambda", "mu")) ->
  p.params
p.params <- p.params[1:5]
as.data.frame(p.params) -> p.params


# set up simulation df

t <- seq(1,100,day_cut)

sim.df <- data.frame(pid = c(1:6),
                     start_state = rep(c(0,1),3),
                     abx_cpt = rep(0,6),
                     tb_start = rep(-999,6),
                     hosp_days = rep(7, 6),
                     abx_days = c(0,0,2,2,7,7))

sim.df$p0 <- as.numeric(sim.df$start_state == 0)
sim.df$p1 <- as.numeric(sim.df$start_state == 1)
sim.df$abx_start <- 0
sim.df$abx_stop <- sim.df$abx_start + sim.df$abx_days
sim.df$abx_stop <- ifelse(sim.df$tb_start == -999,
                          yes = sim.df$abx_stop , no = 1000 )
sim.df$prev_abx <- 999
sim.df$abx_start[sim.df$abx_days == 1] <- -999
sim.df$abx_stop[sim.df$abx_days == 1] <- -999
sim.df$abx_days[sim.df$abx_days == 1] <- 0
sim.df$abx_days[sim.df$abx_days == 1] <- 0
sim.df$abx_stop <- ifelse(sim.df$abx_cpt == 0,yes = sim.df$abx_stop , no = 1000 )
sim.df$hosp_start <- 0
sim.df$hosp_stop <- sim.df$hosp_days
sim.df$prev_hosp <- 999
sim.df$hosp_start[sim.df$hosp_days == 1] <- -999
sim.df$hosp_stop[sim.df$hosp_days == 1] <- -999
sim.df$hosp_days[sim.df$hosp_days == 1] <- 0
sim.df$hosp_days[sim.df$hosp_days == 1] <- 0



# iterate over posterior to solve odes to get states at time t

purrr::pmap(sim.df[,c( "p0", "p1","abx_start", "abx_stop", "prev_abx",
              "hosp_start", "hosp_stop", "prev_hosp","pid")],
     ~iterate_over_posterior_parms(p.params = p.params  ,t = t,
                                   cov_mat = c(..3,..4,..5,..6,..7,..8),
                                   states  =c(..1,..2),
                                   pid = ..9)) -> df.out

do.call(rbind, df.out) -> df.out

# overall estimated prevalance is a weighted mean of those in state one at
# t = 0 and those is state two
# ie for 0.5 at start, p0 + p1 /2
# make this df

left_join(
  df.out %>%
    filter(pid %in% c(1, 3, 5)) %>%
    transmute(
      pid = pid,
      draw = draw,
      time = time,
      pr_esbl_pos_t0esblneg = `2`
    ),
  df.out %>%
    filter(pid %in% c(2, 4, 6)) %>%
    transmute(
      pid = case_when(
        pid == 2 ~ 1,
        pid == 4 ~ 3,
        pid == 6 ~ 5),
      draw = draw,
      time = time,
      pr_esbl_pos_t0esblpos = `2`
    ),
  by = c("pid", "draw", "time")
) %>%
  mutate(pr_esbl_pos =
           (pr_esbl_pos_t0esblneg + pr_esbl_pos_t0esblpos) / 2
  ) %>%
  # link back in the metadata
  left_join(
    sim.df %>%
      select(pid,
             hosp_days,
             abx_days),
    by = "pid"
  ) %>%
  rename(sim_run = pid) ->
  btESBL_model2simulations

# if you want to save
# saveRDS(btESBL_model2simulations, here("data-raw/btESBL_model2simulations.rda"))


