library(rstan)
library(blantyreESBL)
library(deSolve)
library(dplyr)

## simulations of different scenarios from teh posterior of the sitted models

fit_mod2 <- btESBL_model2posterior

p.params <- rstan::extract(fit_mod2)
p.params[1:5] -> p.params
as.data.frame(p.params) -> p.params



# Make stan functions available

expose_stan_functions(file = system.file("extdata",
                                         "ESBLmod_finalV1.0_rk45.stan",
                                         package = "blantyreESBL"),
                      fit_mod2)

# function for deSolve to solve differential eqs

ode_finalstateprob <- function(t, state, parameters) {
  coefs <-
    return_time_varying_coefs_exp_flat (parameters$cov_mat,
                                        t,
                                        parameters$covs_type,
                                        parameters$gammas)
  # print(coefs)

  dp0 <-
    -(parameters$lambda * state[[1]] * exp(sum(parameters$betas * coefs)))   +
    (parameters$mu * state[[2]] * exp(sum(parameters$alphas * coefs)))
  dp1 <-
    (parameters$lambda * state[[1]] * exp(sum(parameters$betas * coefs)))  -
    (parameters$mu * state[[2]] * exp(sum(parameters$alphas * coefs)))
  return(list(c(dp0, dp1)))
}

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

# get the posterior parameters to use

rstan::extract(fit_mod2,
               pars = c("alphas", "betas", "gammas", "lambda", "mu")) ->
  p.params
p.params <- p.params[1:5]
as.data.frame(p.params) -> p.params


# set up simulation df

t <- seq(1,100,day_cut)

sim.df <- data.frame(pid = c(1:3),
                     start_state = rep(1,3),
                     abx_cpt = rep(0,3),
                     tb_start = rep(-999,3),
                     hosp_days = rep(7, 3),
                     abx_days = c(0,2,7))

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
sim.df$p0 <- 0.5
sim.df$p1 <- 0.5


# iterate over posterior to solve odes

purrr::pmap(sim.df[,c( "p0", "p1","abx_start", "abx_stop", "prev_abx",
              "hosp_start", "hosp_stop", "prev_hosp","pid")],
     ~iterate_over_posterior_parms(p.params = p.params  ,t = t,
                                   cov_mat = c(..3,..4,..5,..6,..7,..8),
                                   states  =c(..1,..2),
                                   pid = ..9)) -> df.out

do.call(rbind, df.out) -> df.out

outsum <- df.out %>%
  group_by(time, pid) %>%
  dplyr::summarise(median = median(`2`),
                   lq = quantile(`2`, 0.025)[[1]],
                   uq = quantile(`2`, 0.975)[[1]])

# generate final df

btESBL_model2simulations <-
  left_join (df.out, sim.df, by = "pid")
#write_csv(out.sim.df, "csimulations2.csv")
