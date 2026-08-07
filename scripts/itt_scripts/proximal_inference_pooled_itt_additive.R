##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Proximal Inference Analysis Pooled -- ADDITIVE hazards (exact marginal CIF)
##----- Intention-to-treat
##----- last updated 2026-08-07
##-----
##----- ITT twin of proximal_inference_pooled_pp_additive.R (uses ITT outcomes).


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


# Stage 1: additive NCO (test-negative) model, incl. NCE (flu_vax) --------

prox_add_itt_s1 <- speedglm(Y_neg ~ ns(time_end, knots = c(10,20,30,40,50))*(treatment +
                              sex_admin + age_years + bmi + race + charlson_cat_fac +
                              ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
                              flu_vax),
                            data = dat.long.itt, family = gaussian())
saveRDS(prox_add_itt_s1, paste0(res_path, "prox_add_itt_s1.rds"))

ref_week <- data.frame(time_end = seq(1,53,1), treatment = 0,
                       sex_admin = factor("F"), age_years = 0, bmi = 0,
                       race = factor("White"), charlson_cat_fac = factor("0"),
                       ndi = 0, prior_inf = 0, tests_count = 0,
                       service_region = factor("Central valley"),
                       last_vax_infect_weeks = 0, flu_vax = 0)
ref_week$baseline1 <- predict(prox_add_itt_s1, newdata = ref_week)

dat.long.itt$fitted1 <- predict(prox_add_itt_s1, newdata = dat.long.itt)
dat.long.itt <- left_join(dat.long.itt, ref_week[, c("time_end","baseline1")], by = "time_end")
dat.long.itt$mu <- dat.long.itt$fitted1 - dat.long.itt$baseline1


# Stage 2: additive primary model with bridge; observed primary model ------

prox_add_itt_s2 <- speedglm(Y_pos ~ ns(time_end, knots = c(10,20,30,40,50))*treatment +
                              sex_admin + age_years + bmi + race + charlson_cat_fac +
                              ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
                              mu,
                            data = dat.long.itt, family = gaussian())
saveRDS(prox_add_itt_s2, paste0(res_path, "prox_add_itt_s2.rds"))

prox_add_itt_obs <- speedglm(Y_pos ~ ns(time_end, knots = c(10,20,30,40,50))*treatment +
                               sex_admin + age_years + bmi + race + charlson_cat_fac +
                               ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
                               flu_vax,
                             data = dat.long.itt, family = gaussian())
saveRDS(prox_add_itt_obs, paste0(res_path, "prox_add_itt_obs.rds"))


# De-biased additive causal effect beta_2A(t) -----------------------------

df_ref_A1 <- data.frame(time_end = seq(1,53,1), treatment = 1,
                        sex_admin = factor("F"), age_years = 0, bmi = 0,
                        race = factor("White"), charlson_cat_fac = factor("0"),
                        ndi = 0, prior_inf = 0, tests_count = 0,
                        service_region = factor("Central valley"),
                        last_vax_infect_weeks = 0, mu = 0)
df_ref_A0 <- df_ref_A1; df_ref_A0$treatment <- 0

beta2A <- predict(prox_add_itt_s2, newdata = df_ref_A1) - predict(prox_add_itt_s2, newdata = df_ref_A0)
time_df <- data.frame(time_end = seq(1,53,1), beta2A = beta2A)


# ITT Survival and Risk (exact additive marginal CIF) --------------------

dat$gmaxt <- 53

g <- dat[rep(1:nrow(dat), dat$gmaxt),]
g$time_start <- ave(g$fake_mrn, g$fake_mrn, FUN=seq_along)
g$time_start <- (g$time_start-1)*time_unit
g$time_end   <- g$time_start + time_unit
g$treatment_obs <- g$treatment

g$lambda1 <- predict(prox_add_itt_s1,  newdata = g)
g$lambda2 <- predict(prox_add_itt_obs, newdata = g)
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

pci.add.itt.risk.pointest <- tibble(
  sim = 0L,
  time_end = prox.add.res$time_end,
  risk0 = prox.add.res$risk0,
  risk1 = prox.add.res$risk1
)

saveRDS(pci.add.itt.risk.pointest, paste0(res_path, "pci.add.itt.risk.pointest.rds"))
