## simulations of different scenarios
# plan: n_pat patients
#prop cpt: prop_cpt {0.1-0.5}
# prop tb: prop_tb {0.1-0.5] started at a rate 1/tb_day per day such that prop_tb pts are on tb rx by
# day prop_tb
# hosp_adm_days {2,5,7,10,14} days
# abx_days {2,5,7,10,14} days

# day_cut is no of days between prevalence estimate
# start prevalence ; start_esbl_prop = 0.5

# plan - generate a df with one row for each participant with covariate vals for that pt
# then expand to tstart- tstop and build up states using ode solver

fit_mod2 <- readRDS("chapter_9/stan_models/model_2/stanfit_m2.rds")

p.params <- rstan::extract(fit_mod2)
p.params[1:5] -> p.params
as.data.frame(p.params) -> p.params

library(tidyverse)
library(rstan)
library(deSolve)
expose_stan_functions("other_scripts/stan_model/ESBLmod_finalV1.0_rk45.stan")


ode_finalstateprob <- function(t, state, parameters) {
  coefs <- return_time_varying_coefs_exp_flat (parameters$cov_mat, t, parameters$covs_type, parameters$gammas)
  # print(coefs)
  
  dp0 <- -(parameters$lambda * state[[1]] * exp(sum(parameters$betas * coefs)))   + (parameters$mu * state[[2]] * exp(sum(parameters$alphas * coefs)))
  dp1 <- (parameters$lambda * state[[1]] * exp(sum(parameters$betas *coefs)))  - (parameters$mu * state[[2]] * exp(sum(parameters$alphas * coefs)))
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


rstan::extract(fit_mod2, pars = c("alphas", "betas", "gammas", "lambda", "mu")) -> p.params
p.params <- p.params[1:5]
as.data.frame(p.params) -> p.params

t <- seq(1,100,day_cut)

sim.df <- data.frame(pid = c(1:3), start_state = rep(1,3), abx_cpt = rep(0,3), tb_start = rep(-999,3), hosp_days = rep(7, 3), abx_days = c(0,2,7))






sim.df$p0 <- as.numeric(sim.df$start_state ==0)
sim.df$p1 <- as.numeric(sim.df$start_state ==1)
sim.df$abx_start <- 0
sim.df$abx_stop <- sim.df$abx_start + sim.df$abx_days
sim.df$abx_stop <- ifelse(sim.df$tb_start == -999,yes = sim.df$abx_stop , no = 1000 )

sim.df$prev_abx <- 999
sim.df$abx_start[sim.df$abx_days == 1] <- -999
sim.df$abx_stop[sim.df$abx_days == 1] <- -999
sim.df$abx_days[sim.df$abx_days == 1] <- 0
sim.df$abx_days[sim.df$abx_days == 1] <- 0
sim.df$abx_stop <- ifelse(sim.df$abx_cpt== 0,yes = sim.df$abx_stop , no = 1000 )

sim.df$hosp_start <- 0
sim.df$hosp_stop <- sim.df$hosp_days
sim.df$prev_hosp <- 999
sim.df$hosp_start[sim.df$hosp_days == 1] <- -999
sim.df$hosp_stop[sim.df$hosp_days == 1] <- -999
sim.df$hosp_days[sim.df$hosp_days == 1] <- 0
sim.df$hosp_days[sim.df$hosp_days == 1] <- 0


sim.df$p0 <- 0.5
sim.df$p1 <- 0.5



purrr::pmap(sim.df[,c( "p0", "p1","abx_start", "abx_stop", "prev_abx",
              "hosp_start", "hosp_stop", "prev_hosp","pid")],
     ~iterate_over_posterior_parms(p.params = p.params  ,t = t,
                                   cov_mat = c(..3,..4,..5,..6,..7,..8), states  =c(..1,..2),
                                   pid = ..9)) -> test

do.call(rbind, test) -> test
test
outsum <- test %>% group_by(time, pid) %>% dplyr::summarise(median = median(`2`),
                                               lq = quantile(`2`, 0.025)[[1]],
                                               uq = quantile(`2`, 0.975)[[1]])

out.sim.df <- left_join (test, sim.df, by = "pid")
write_csv(out.sim.df, "chapter_9/simulations2.csv")

outsum <- merge(outsum, select(sim.df, pid, hosp_days, abx_days, abx_cpt), all.x =T)

outsum$hosp_days_str <- paste0("Hosp: ",outsum$hosp_days, "d" )
outsum$hosp_days_str <- factor(outsum$hosp_days_str, levels = unique(outsum$hosp_days_str[order(outsum$hosp_days)]))
outsum$abx_days_str <- paste0(outsum$abx_days, " days" )
outsum$abx_days_str <- factor(outsum$abx_days_str, levels = unique(outsum$abx_days_str[order(outsum$abx_days)]))

outsum$abx_cpt[outsum$abx_cpt == 0] <- "No CPT"
outsum$abx_cpt[outsum$abx_cpt == 1] <- "CPT"

#write.csv(outsum, "chapter_9/simulations.csv", row.names = FALSE)

ggplot(outsum, aes(time, median, ymin = lq, ymax = uq, group = abx_cpt, linetype = as.factor(abx_cpt))) +
  geom_line() + geom_ribbon(alpha = 0.2, color = NA) + facet_grid(hosp_days_str ~ abx_days_str) + theme_bw() + theme(legend.title = element_blank(), panel.spacing = unit(1, "line")) + scale_linetype_manual(values = c("dotted", "solid")) + xlab("Time (days") + ylab("Pr(ESBL)")

# color = as.factor(abx_cpt), fill = as.factor(abx_cpt),

outsum$abx_days_str <- as.factor(outsum$abx_days_str)

ggplot(outsum, aes(time, median, ymin = lq, ymax = uq,group = fct_rev(abx_days_str), color = fct_rev(abx_days_str), fill = fct_rev(abx_days_str))) +
  geom_line() + geom_ribbon(alpha = 0.2, color = NA) +
  #facet_grid(hosp_days_str ~ abx_days_str) + 
  theme_bw() + theme(text = element_text(size =14)) + xlab("Time (days)") + ylab("Pr(ESBL)") +
  scale_color_manual(name = "Days of antibiotic", values = c("#00468BFF", "#ED0000FF", "#42B540FF")) +
 scale_fill_manual(name = "Days of antibiotic",values = c("#00468BFF", "#ED0000FF", "#42B540FF"))
  