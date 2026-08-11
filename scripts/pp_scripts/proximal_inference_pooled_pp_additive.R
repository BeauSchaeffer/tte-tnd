##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Proximal Inference Analysis Pooled -- ADDITIVE hazards (exact marginal CIF)
##----- Per-protocol, no censoring weights
##----- last updated 2026-08-07
##-----
##----- Additive analog of proximal_inference_pooled_pp.R (Li et al., arXiv:2409.08924,
##----- "competing risks as negative controls"). Pooled LINEAR (identity-link)
##----- cause-specific hazard models replace the pooled logistic models. The bridge
##----- is the fitted stage-1 hazard on the identity scale (not its log), and the
##----- intervention enters as an additive shift beta_2A(t)*(a - A), which factors
##----- out of the exp(-cumulative hazard) survival -- a frailty-exact marginal CIF.


# Packages ----------------------------------------------------------------


library(tidyverse)
library(data.table)
library(speedglm)
library(splines)


# Data --------------------------------------------------------------------


data_Y3 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch.3/data_Y3_weekmatch.rds")
dat <- data_Y3
setDT(dat)

res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results_pp.5/"


# PP long format expansion -----------------------------------------------


time_unit <- 1

dat$max_units <- ceiling(dat$Y3_pp_t_trunc/time_unit)+1
dat.long.pp <- dat[rep(1:nrow(dat), dat$max_units),]

dat.long.pp$time_start <- ave(dat.long.pp$fake_mrn, dat.long.pp$fake_mrn, FUN=seq_along)
dat.long.pp$time_start <- (dat.long.pp$time_start-1)*time_unit
dat.long.pp$time_end <- dat.long.pp$time_start+time_unit

dat.long.pp$Y_pos <- ifelse(dat.long.pp$Y3_pp_trunc == 2 & dat.long.pp$Y3_pp_t_trunc == dat.long.pp$time_start, 1, 0)
dat.long.pp$Y_neg <- ifelse(dat.long.pp$Y3_pp_trunc == 1 & dat.long.pp$Y3_pp_t_trunc == dat.long.pp$time_start, 1, 0)
dat.long.pp$C     <- ifelse(dat.long.pp$Y3_pp_trunc == 0 & dat.long.pp$Y3_pp_t_trunc == dat.long.pp$time_start, 1, 0)

dat.long.pp$Y_pos <- ifelse(dat.long.pp$C==1, NA, dat.long.pp$Y_pos)
dat.long.pp$Y_neg <- ifelse(dat.long.pp$C==1, NA, dat.long.pp$Y_neg)


# Stage 1: additive NCO (test-negative) model, incl. NCE (flu_vax) ---------

prox_add_pp_s1 <- speedglm(Y_neg ~ ns(time_end, knots = c(10,20,30,40,50))*(treatment +
                             sex_admin + age_years + bmi + race + charlson_cat_fac +
                             ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
                             flu_vax),
                           data = dat.long.pp, family = gaussian(), sparse=FALSE)
saveRDS(prox_add_pp_s1, paste0(res_path, "prox_add_pp_s1.rds"))

## Bridge mu = covariate-driven part of the fitted stage-1 hazard (identity scale).
## We subtract the time-varying baseline (fitted at reference covariates) so mu is
## not collinear with the stage-2 ns(time) baseline. The baseline depends only on
## time, so we evaluate it on a 53-week reference frame and join by time.
ref_week <- data.frame(time_end = seq(1,53,1), treatment = 0,
                       sex_admin = factor("F"), age_years = 0, bmi = 0,
                       race = factor("White"), charlson_cat_fac = factor("0"),
                       ndi = 0, prior_inf = 0, tests_count = 0,
                       service_region = factor("Central valley"),
                       last_vax_infect_weeks = 0, flu_vax = 0)
ref_week$baseline1 <- predict(prox_add_pp_s1, newdata = ref_week)

dat.long.pp$fitted1 <- predict(prox_add_pp_s1, newdata = dat.long.pp)
dat.long.pp <- left_join(dat.long.pp, ref_week[, c("time_end","baseline1")], by = "time_end")
dat.long.pp$mu <- dat.long.pp$fitted1 - dat.long.pp$baseline1


# Stage 2: additive primary (test-positive) model with the bridge ---------
## no NCE (flu_vax) in stage 2; the treatment coefficient is the causal beta_2A(t)

prox_add_pp_s2 <- speedglm(Y_pos ~ ns(time_end, knots = c(10,20,30,40,50))*treatment +
                             sex_admin + age_years + bmi + race + charlson_cat_fac +
                             ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
                             mu,
                           data = dat.long.pp, family = gaussian(), sparse=FALSE)
saveRDS(prox_add_pp_s2, paste0(res_path, "prox_add_pp_s2.rds"))


# Observed primary hazard model (with NCE, no bridge) ---------------------
## used as the anchor for the g-computation, evaluated at observed treatment

prox_add_pp_obs <- speedglm(Y_pos ~ ns(time_end, knots = c(10,20,30,40,50))*treatment +
                              sex_admin + age_years + bmi + race + charlson_cat_fac +
                              ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
                              flu_vax,
                            data = dat.long.pp, family = gaussian(), sparse=FALSE)
saveRDS(prox_add_pp_obs, paste0(res_path, "prox_add_pp_obs.rds"))


# De-biased additive causal effect beta_2A(t) -----------------------------
## stage-2 treatment contrast at reference covariates with the bridge fixed at 0

df_ref_A1 <- data.frame(time_end = seq(1,53,1), treatment = 1,
                        sex_admin = factor("F"), age_years = 0, bmi = 0,
                        race = factor("White"), charlson_cat_fac = factor("0"),
                        ndi = 0, prior_inf = 0, tests_count = 0,
                        service_region = factor("Central valley"),
                        last_vax_infect_weeks = 0, mu = 0)
df_ref_A0 <- df_ref_A1; df_ref_A0$treatment <- 0

beta2A <- predict(prox_add_pp_s2, newdata = df_ref_A1) - predict(prox_add_pp_s2, newdata = df_ref_A0)
time_df <- data.frame(time_end = seq(1,53,1), beta2A = beta2A)


# PP Survival and Risk (exact additive marginal CIF) ----------------------

dat$gmaxt <- 53

g <- dat[rep(1:nrow(dat), dat$gmaxt),]
g$time_start <- ave(g$fake_mrn, g$fake_mrn, FUN=seq_along)
g$time_start <- (g$time_start-1)*time_unit
g$time_end   <- g$time_start + time_unit
g$treatment_obs <- g$treatment

## competing (NCO) hazard from stage 1; observed primary hazard from the observed
## model -- both at the OBSERVED treatment. The intervention enters only via beta2A.
g$lambda1 <- predict(prox_add_pp_s1,  newdata = g)   # test-negative (competing)
g$lambda2 <- predict(prox_add_pp_obs, newdata = g)   # test-positive (primary), observed

g$lambda1 <- pmin(pmax(g$lambda1, 0), 1)
g$lambda2 <- pmin(pmax(g$lambda2, 0), 1)

g <- left_join(g, time_df, by = "time_end")

g <- g |>
  arrange(fake_mrn, time_end) |>
  group_by(fake_mrn) |>
  mutate(
    hz0 = pmin(pmax(lambda2 + beta2A * (0 - treatment_obs), 0), 1),
    hz1 = pmin(pmax(lambda2 + beta2A * (1 - treatment_obs), 0), 1),
    surv0 = exp(-cumsum(pmin(pmax(lambda1 + hz0, 0), 1))),
    surv1 = exp(-cumsum(pmin(pmax(lambda1 + hz1, 0), 1))),
    risk0 = cumsum(hz0 * lag(surv0, default = 1)),
    risk1 = cumsum(hz1 * lag(surv1, default = 1))
  ) |>
  ungroup()

prox.add.res <- g |>
  group_by(time_end) |>
  summarise(risk0 = mean(risk0), risk1 = mean(risk1), .groups = "drop")

pci.add.pp.risk.pointest <- tibble(
  sim = 0L,
  time_end = prox.add.res$time_end,
  risk0 = prox.add.res$risk0,   # marginal CIF under no booster
  risk1 = prox.add.res$risk1    # marginal CIF under booster
)

saveRDS(pci.add.pp.risk.pointest, paste0(res_path, "pci.add.pp.risk.pointest.rds"))
