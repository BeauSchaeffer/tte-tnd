##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Equi-Confounding Analysis Pooled -- ADDITIVE hazards (exact marginal CIF)
##----- Intention-to-treat
##----- last updated 2026-08-07
##-----
##----- ITT twin of equi_confounding_pooled_pp_additive.R (uses ITT outcomes).


# Packages ----------------------------------------------------------------


library(tidyverse)
library(data.table)
library(speedglm)
library(splines)


# Data --------------------------------------------------------------------


data_Y3 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch.3/data_Y3_weekmatch.rds")
dat <- data_Y3
setDT(dat)

res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results_itt.5/"


# ITT long format expansion ----------------------------------------------


time_unit <- 1

dat$max_units <- ceiling(dat$Y3_itt_t_trunc/time_unit)+1
dat.long.itt <- dat[rep(1:nrow(dat), dat$max_units),]

dat.long.itt$time_start <- ave(dat.long.itt$fake_mrn, dat.long.itt$fake_mrn, FUN=seq_along)
dat.long.itt$time_start <- (dat.long.itt$time_start-1)*time_unit
dat.long.itt$time_end <- dat.long.itt$time_start+time_unit

dat.long.itt$Y_pos <- ifelse(dat.long.itt$Y3_itt_trunc == 2 & dat.long.itt$Y3_itt_t_trunc == dat.long.itt$time_start, 1, 0)
dat.long.itt$Y_neg <- ifelse(dat.long.itt$Y3_itt_trunc == 1 & dat.long.itt$Y3_itt_t_trunc == dat.long.itt$time_start, 1, 0)
dat.long.itt$C     <- ifelse(dat.long.itt$Y3_itt_trunc == 0 & dat.long.itt$Y3_itt_t_trunc == dat.long.itt$time_start, 1, 0)

dat.long.itt$Y_pos <- ifelse(dat.long.itt$C==1, NA, dat.long.itt$Y_pos)
dat.long.itt$Y_neg <- ifelse(dat.long.itt$C==1, NA, dat.long.itt$Y_neg)


# ITT Pooled ADDITIVE hazards --------------------------------------------


eqc_add_itt_fit1 <- speedglm(Y_neg ~ ns(time_end, knots = c(10,20,30,40,50))*treatment +
                               sex_admin + age_years + bmi + race + charlson_cat_fac +
                               ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
                               flu_vax,
                             data = dat.long.itt, family = gaussian())

eqc_add_itt_fit2 <- speedglm(Y_pos ~ ns(time_end, knots = c(10,20,30,40,50))*treatment +
                               sex_admin + age_years + bmi + race + charlson_cat_fac +
                               ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
                               flu_vax,
                             data = dat.long.itt, family = gaussian())

saveRDS(eqc_add_itt_fit1, paste0(res_path, "eqc_add_itt_fit1.rds"))
saveRDS(eqc_add_itt_fit2, paste0(res_path, "eqc_add_itt_fit2.rds"))


# De-biased additive causal effect beta_2A(t) -----------------------------


ref_A1 <- data.frame(time_end = seq(1,53,1), treatment = 1,
                     sex_admin = factor("F"), age_years = 0, bmi = 0,
                     race = factor("White"), charlson_cat_fac = factor("0"),
                     ndi = 0, prior_inf = 0, tests_count = 0,
                     service_region = factor("Central valley"),
                     last_vax_infect_weeks = 0, flu_vax = 0)
ref_A0 <- ref_A1; ref_A0$treatment <- 0

delta1 <- predict(eqc_add_itt_fit1, newdata = ref_A1) - predict(eqc_add_itt_fit1, newdata = ref_A0)
delta2 <- predict(eqc_add_itt_fit2, newdata = ref_A1) - predict(eqc_add_itt_fit2, newdata = ref_A0)

time_df <- data.frame(time_end = seq(1,53,1), beta2A = delta2 - delta1)


# ITT Survival and Risk (exact additive marginal CIF) --------------------


dat$gmaxt <- 53

g <- dat[rep(1:nrow(dat), dat$gmaxt),]
g$time_start <- ave(g$fake_mrn, g$fake_mrn, FUN=seq_along)
g$time_start <- (g$time_start-1)*time_unit
g$time_end   <- g$time_start + time_unit
g$treatment_obs <- g$treatment

g$lambda1 <- predict(eqc_add_itt_fit1, newdata = g)
g$lambda2 <- predict(eqc_add_itt_fit2, newdata = g)
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

eqc.add.res <- g |>
  group_by(time_end) |>
  summarise(risk0 = mean(risk0), risk1 = mean(risk1), .groups = "drop")

eqc.add.itt.risk.pointest <- tibble(
  sim = 0L,
  time_end = eqc.add.res$time_end,
  risk0 = eqc.add.res$risk0,
  risk1 = eqc.add.res$risk1
)

saveRDS(eqc.add.itt.risk.pointest, paste0(res_path, "eqc.add.itt.risk.pointest.rds"))
