##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Equi-Confounding Analysis Pooled -- ADDITIVE hazards (exact marginal CIF)
##----- Per-protocol, no censoring weights
##----- last updated 2026-08-07
##-----
##----- Additive analog of equi_confounding_pooled_pp.R.
##----- Pooled LINEAR (identity-link) cause-specific hazard models replace the
##----- pooled logistic models: the fitted mean IS the additive discrete-time
##----- hazard. The intervention enters as an additive shift beta_2A(t)*(a - A),
##----- which is U-independent and therefore factors out of the survival built as
##----- exp(-cumulative hazard) -- giving a frailty-exact marginal CIF.


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


# PP Pooled ADDITIVE hazards ----------------------------------------------

## family = gaussian() => identity link: the fitted mean is the additive
## discrete-time cause-specific hazard (a linear probability model per person-week).

eqc_add_pp_fit1 <- speedglm(Y_neg ~ ns(time_end, knots = c(10,20,30,40,50))*treatment +
                              sex_admin + age_years + bmi + race + charlson_cat_fac +
                              ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
                              flu_vax,
                            data = dat.long.pp, family = gaussian())

eqc_add_pp_fit2 <- speedglm(Y_pos ~ ns(time_end, knots = c(10,20,30,40,50))*treatment +
                              sex_admin + age_years + bmi + race + charlson_cat_fac +
                              ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
                              flu_vax,
                            data = dat.long.pp, family = gaussian())

saveRDS(eqc_add_pp_fit1, paste0(res_path, "eqc_add_pp_fit1.rds"))
saveRDS(eqc_add_pp_fit2, paste0(res_path, "eqc_add_pp_fit2.rds"))


# De-biased additive causal effect beta_2A(t) -----------------------------

## Equi-confounding subtraction on the identity (hazard-difference) scale:
##   beta_2A(t) = [lambda2(A=1)-lambda2(A=0)] - [lambda1(A=1)-lambda1(A=0)]
## evaluated at reference covariates as a function of time. Non-treatment
## covariates cancel in each within-model difference.

ref_A1 <- data.frame(time_end = seq(1,53,1), treatment = 1,
                     sex_admin = factor("F"), age_years = 0, bmi = 0,
                     race = factor("White"), charlson_cat_fac = factor("0"),
                     ndi = 0, prior_inf = 0, tests_count = 0,
                     service_region = factor("Central valley"),
                     last_vax_infect_weeks = 0, flu_vax = 0)
ref_A0 <- ref_A1; ref_A0$treatment <- 0

delta1 <- predict(eqc_add_pp_fit1, newdata = ref_A1) - predict(eqc_add_pp_fit1, newdata = ref_A0)  # NCO (confounding) diff
delta2 <- predict(eqc_add_pp_fit2, newdata = ref_A1) - predict(eqc_add_pp_fit2, newdata = ref_A0)  # observed test-pos diff

time_df <- data.frame(time_end = seq(1,53,1), beta2A = delta2 - delta1)


# PP Survival and Risk (exact additive marginal CIF) ----------------------

dat$gmaxt <- 53

## single g-formula frame; observed cause-specific hazards are evaluated at the
## OBSERVED treatment, and the intervention enters only through the additive switch.
g <- dat[rep(1:nrow(dat), dat$gmaxt),]
g$time_start <- ave(g$fake_mrn, g$fake_mrn, FUN=seq_along)
g$time_start <- (g$time_start-1)*time_unit
g$time_end   <- g$time_start + time_unit
g$treatment_obs <- g$treatment

g$lambda1 <- predict(eqc_add_pp_fit1, newdata = g)   # test-negative (competing / NCO)
g$lambda2 <- predict(eqc_add_pp_fit2, newdata = g)   # test-positive (primary), observed

## additive linear-probability hazards can fall outside [0,1]; clamp
g$lambda1 <- pmin(pmax(g$lambda1, 0), 1)
g$lambda2 <- pmin(pmax(g$lambda2, 0), 1)

g <- left_join(g, time_df, by = "time_end")

## counterfactual curves under a=0 and a=1 via additive switch + exp-form survival
g <- g |>
  arrange(fake_mrn, time_end) |>
  group_by(fake_mrn) |>
  mutate(
    ## additive switch on the primary (test-positive) hazard: lambda2 + beta2A*(a - A)
    hz0 = pmin(pmax(lambda2 + beta2A * (0 - treatment_obs), 0), 1),
    hz1 = pmin(pmax(lambda2 + beta2A * (1 - treatment_obs), 0), 1),
    ## all-cause hazard (competing risk lambda1 is unaffected by treatment)
    ## exp-form survival makes the additive shift factor out (frailty-exact)
    surv0 = exp(-cumsum(pmin(pmax(lambda1 + hz0, 0), 1))),
    surv1 = exp(-cumsum(pmin(pmax(lambda1 + hz1, 0), 1))),
    ## Aalen-Johansen sub-density: hazard * survival to start of the interval
    risk0 = cumsum(hz0 * lag(surv0, default = 1)),
    risk1 = cumsum(hz1 * lag(surv1, default = 1))
  ) |>
  ungroup()

## standardize (average) over the observed covariate/treatment distribution
eqc.add.res <- g |>
  group_by(time_end) |>
  summarise(risk0 = mean(risk0), risk1 = mean(risk1), .groups = "drop")

eqc.add.pp.risk.pointest <- tibble(
  sim = 0L,
  time_end = eqc.add.res$time_end,
  risk0 = eqc.add.res$risk0,   # marginal CIF under no booster
  risk1 = eqc.add.res$risk1    # marginal CIF under booster
)

saveRDS(eqc.add.pp.risk.pointest, paste0(res_path, "eqc.add.pp.risk.pointest.rds"))
