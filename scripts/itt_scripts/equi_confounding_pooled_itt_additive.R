##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Equi-Confounding Analysis Pooled -- ADDITIVE hazards
##----- Intention-to-treat
##----- last updated 2026-08-12
##-----
##----- ITT twin of equi_confounding_pooled_pp_additive.R (uses ITT outcomes),
##----- kept line-for-line parallel to equi_confounding_pooled_itt.R.


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


# ITT Survival and Risk (parallel to multiplicative; ADDITIVE de-biasing) -


dat$gmaxt <- 53

### G formula data setup A=0
eqc_itt_A0.long <- dat[rep(1:nrow(dat), dat$gmaxt),]
eqc_itt_A0.long$time_start <- ave(eqc_itt_A0.long$fake_mrn, eqc_itt_A0.long$fake_mrn, FUN=seq_along)
eqc_itt_A0.long$time_start <- (eqc_itt_A0.long$time_start-1)*time_unit
eqc_itt_A0.long$time_end <- eqc_itt_A0.long$time_start+time_unit
eqc_itt_A0.long$treatment <- 0

### G formula data setup A=1
eqc_itt_A1.long <- dat[rep(1:nrow(dat), dat$gmaxt),]
eqc_itt_A1.long$time_start <- ave(eqc_itt_A1.long$fake_mrn, eqc_itt_A1.long$fake_mrn, FUN=seq_along)
eqc_itt_A1.long$time_start <- (eqc_itt_A1.long$time_start-1)*time_unit
eqc_itt_A1.long$time_end <- eqc_itt_A1.long$time_start+time_unit
eqc_itt_A1.long$treatment <- 1

### Calculate predicted hazards (identity link -> additive discrete-time hazard):
eqc_itt_A0.long$hazard_pos <- predict(eqc_add_itt_fit2, newdata=eqc_itt_A0.long)
eqc_itt_A1.long$hazard_pos <- predict(eqc_add_itt_fit2, newdata=eqc_itt_A1.long)
eqc_itt_A0.long$hazard_neg <- predict(eqc_add_itt_fit1, newdata=eqc_itt_A0.long)
eqc_itt_A1.long$hazard_neg <- predict(eqc_add_itt_fit1, newdata=eqc_itt_A1.long)
### identity-link hazards are not range-restricted; clamp to [0,1]
eqc_itt_A0.long$hazard_pos <- pmin(pmax(eqc_itt_A0.long$hazard_pos, 0), 1)
eqc_itt_A1.long$hazard_pos <- pmin(pmax(eqc_itt_A1.long$hazard_pos, 0), 1)
eqc_itt_A0.long$hazard_neg <- pmin(pmax(eqc_itt_A0.long$hazard_neg, 0), 1)
eqc_itt_A1.long$hazard_neg <- pmin(pmax(eqc_itt_A1.long$hazard_neg, 0), 1)
### Corrected hazards under no treatment -- ADDITIVE de-biasing (cf. multiplicative ratio):
### lambda2(A=0) + [lambda1(A=1) - lambda1(A=0)]
eqc_itt_A0.long$hazard_pos_c <- eqc_itt_A0.long$hazard_pos + (eqc_itt_A1.long$hazard_neg - eqc_itt_A0.long$hazard_neg)
eqc_itt_A0.long$hazard_pos_c <- pmin(pmax(eqc_itt_A0.long$hazard_pos_c, 0), 1)

### Calculate (1 - hazard)
eqc_itt_A0.long$pnoevent_pos <- 1 - eqc_itt_A0.long$hazard_pos
eqc_itt_A1.long$pnoevent_pos <- 1 - eqc_itt_A1.long$hazard_pos
eqc_itt_A0.long$pnoevent_neg <- 1 - eqc_itt_A0.long$hazard_neg
eqc_itt_A1.long$pnoevent_neg <- 1 - eqc_itt_A1.long$hazard_neg
### Corrected (1 - hazard) under no treatment
eqc_itt_A0.long$pnoevent_pos_c <- 1 - eqc_itt_A0.long$hazard_pos_c
### corrected curve uses the treated population's test-negative (competing) hazard,
### invariant to treatment under the NCO assumption (lambda_1 at A=1)
eqc_itt_A0.long$pnoevent_neg_c <- 1 - eqc_itt_A1.long$hazard_neg

### Sort the data by ID, time
eqc_itt_A0.long <- eqc_itt_A0.long[order(eqc_itt_A0.long$fake_mrn, eqc_itt_A0.long$time_end),]
eqc_itt_A1.long <- eqc_itt_A1.long[order(eqc_itt_A1.long$fake_mrn, eqc_itt_A1.long$time_end),]

### Calculate the cumulative survival

# lag P(no event pos)

eqc_itt_A0.long <- eqc_itt_A0.long |>
  arrange(fake_mrn, time_end) |>
  group_by(fake_mrn) |>
  mutate(pnoevent_pos_lag = lag(pnoevent_pos, n=1, default=1),
         pnoevent_pos_c_lag = lag(pnoevent_pos_c, n=1, default=1)) |>
  ungroup()

eqc_itt_A1.long <- eqc_itt_A1.long |>
  arrange(fake_mrn, time_end) |>
  group_by(fake_mrn) |>
  mutate(pnoevent_pos_lag = lag(pnoevent_pos, n=1, default=1)) |>
  ungroup()

# product at each time (P(no event neg) * lag P(no event pos))

eqc_itt_A0.long$surv_prod_lag <- eqc_itt_A0.long$pnoevent_neg * eqc_itt_A0.long$pnoevent_pos_lag
eqc_itt_A1.long$surv_prod_lag <- eqc_itt_A1.long$pnoevent_neg * eqc_itt_A1.long$pnoevent_pos_lag
eqc_itt_A0.long$surv_prod_c_lag <- eqc_itt_A0.long$pnoevent_neg_c * eqc_itt_A0.long$pnoevent_pos_c_lag

# cumulative product within individual

eqc_itt_A0.long$survival_pos <- ave(eqc_itt_A0.long$surv_prod_lag, eqc_itt_A0.long$fake_mrn, FUN=cumprod)
eqc_itt_A1.long$survival_pos <- ave(eqc_itt_A1.long$surv_prod_lag, eqc_itt_A1.long$fake_mrn, FUN=cumprod)
eqc_itt_A0.long$survival_pos_c <- ave(eqc_itt_A0.long$surv_prod_c_lag, eqc_itt_A0.long$fake_mrn, FUN=cumprod)

### Calculate risk using CIF estimator

# product at each time (haz pos * surv pos)

eqc_itt_A0.long$risk_prod_pos <- eqc_itt_A0.long$hazard_pos * eqc_itt_A0.long$survival_pos
eqc_itt_A1.long$risk_prod_pos <- eqc_itt_A1.long$hazard_pos * eqc_itt_A1.long$survival_pos
eqc_itt_A0.long$risk_prod_pos_c <- eqc_itt_A0.long$hazard_pos_c * eqc_itt_A0.long$survival_pos_c

# cumulative sum within individual

eqc_itt_A0.long$risk_pos <- ave(eqc_itt_A0.long$risk_prod_pos, eqc_itt_A0.long$fake_mrn, FUN=cumsum)
eqc_itt_A1.long$risk_pos <- ave(eqc_itt_A1.long$risk_prod_pos, eqc_itt_A1.long$fake_mrn, FUN=cumsum)
eqc_itt_A0.long$risk_pos_c <- ave(eqc_itt_A0.long$risk_prod_pos_c, eqc_itt_A0.long$fake_mrn, FUN=cumsum)

# Calculate the average risk at each time point

eqc_itt_A0.long.res <- aggregate(risk_pos ~ time_end, data=eqc_itt_A0.long, FUN=mean)
eqc_itt_A1.long.res <- aggregate(risk_pos ~ time_end, data=eqc_itt_A1.long, FUN=mean)
eqc_itt_A0.long.res.c <- aggregate(risk_pos_c ~ time_end, data=eqc_itt_A0.long, FUN=mean)

# save point estimate risk curves in bootstrap-compatible format
eqc.add.itt.risk.pointest <- tibble(
  sim = 0L,  # 0 = main analysis (bootstraps are 1..B)
  time_end = eqc_itt_A0.long.res$time_end,
  risk0 = eqc_itt_A0.long.res$risk_pos,
  risk0corr = eqc_itt_A0.long.res.c$risk_pos_c,
  risk1 = eqc_itt_A1.long.res$risk_pos
)

saveRDS(eqc.add.itt.risk.pointest, paste0(res_path, "eqc.add.itt.risk.pointest.rds"))
