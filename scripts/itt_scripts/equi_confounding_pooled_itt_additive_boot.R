##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Equi-Confounding Pooled ADDITIVE ITT Bootstrap
##----- Intention-to-treat
##----- last updated 2026-08-12
##-----
##----- ITT twin of equi_confounding_pooled_pp_additive_boot.R, kept line-for-line
##----- parallel to equi_confounding_pooled_itt.R / _boot.R.


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


# Boot setup --------------------------------------------------------------


num.boot <- 200

set.seed(1155)
seed <- floor(runif(num.boot)*10^8)

setkey(dat, subclass)
subclasses <- dat[, unique(subclass)]
n_sub <- length(subclasses)

time_unit <- 1


boot.results <- lapply(1:num.boot, function(i){

  t0 <- Sys.time()
  set.seed(seed[i])
  message("Starting EQC ITT additive bootstrap ", i, " (seed=", seed[i], ")")

  ## resample matched sets, unique ids per draw
  samp_sub <- sample(subclasses, size = n_sub, replace = TRUE)
  map <- data.table(j = seq_along(samp_sub), subclass = samp_sub)
  dat.boot <- dat[map, on = "subclass", allow.cartesian = TRUE]
  dat.boot[, bootid := j]
  dat.boot[, bootid_mrn := .GRP, by = .(bootid, fake_mrn)]

  ## long format (person-week), keyed on bootid_mrn
  dat.boot$max_units <- ceiling(dat.boot$Y3_itt_t_trunc/time_unit)+1
  dat.long <- dat.boot[rep(1:nrow(dat.boot), dat.boot$max_units),]
  dat.long$time_start <- ave(dat.long$bootid_mrn, dat.long$bootid_mrn, FUN=seq_along)
  dat.long$time_start <- (dat.long$time_start-1)*time_unit
  dat.long$time_end   <- dat.long$time_start+time_unit
  dat.long <- dat.long[order(dat.long$bootid_mrn, dat.long$time_end),]
  dat.long$Y_pos <- ifelse(dat.long$Y3_itt_trunc==2 & dat.long$Y3_itt_t_trunc==dat.long$time_start, 1, 0)
  dat.long$Y_neg <- ifelse(dat.long$Y3_itt_trunc==1 & dat.long$Y3_itt_t_trunc==dat.long$time_start, 1, 0)
  dat.long$C     <- ifelse(dat.long$Y3_itt_trunc==0 & dat.long$Y3_itt_t_trunc==dat.long$time_start, 1, 0)
  dat.long$Y_pos <- ifelse(dat.long$C==1, NA, dat.long$Y_pos)
  dat.long$Y_neg <- ifelse(dat.long$C==1, NA, dat.long$Y_neg)

  ## additive (identity-link) cause-specific hazard models
  fit1 <- speedglm(Y_neg ~ ns(time_end, knots = c(10,20,30,40,50))*treatment +
                     sex_admin + age_years + bmi + race + charlson_cat_fac +
                     ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
                     flu_vax,
                   data = dat.long, family = gaussian())
  fit2 <- speedglm(Y_pos ~ ns(time_end, knots = c(10,20,30,40,50))*treatment +
                     sex_admin + age_years + bmi + race + charlson_cat_fac +
                     ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
                     flu_vax,
                   data = dat.long, family = gaussian())

  ## g formula setup
  dat.boot$gmaxt <- 53

  eqc_itt_A0.long <- dat.boot[rep(1:nrow(dat.boot), dat.boot$gmaxt),]
  eqc_itt_A0.long$time_start <- ave(eqc_itt_A0.long$bootid_mrn, eqc_itt_A0.long$bootid_mrn, FUN=seq_along)
  eqc_itt_A0.long$time_start <- (eqc_itt_A0.long$time_start-1)*time_unit
  eqc_itt_A0.long$time_end <- eqc_itt_A0.long$time_start+time_unit
  eqc_itt_A0.long$treatment <- 0

  eqc_itt_A1.long <- dat.boot[rep(1:nrow(dat.boot), dat.boot$gmaxt),]
  eqc_itt_A1.long$time_start <- ave(eqc_itt_A1.long$bootid_mrn, eqc_itt_A1.long$bootid_mrn, FUN=seq_along)
  eqc_itt_A1.long$time_start <- (eqc_itt_A1.long$time_start-1)*time_unit
  eqc_itt_A1.long$time_end <- eqc_itt_A1.long$time_start+time_unit
  eqc_itt_A1.long$treatment <- 1

  ### Calculate predicted hazards (identity link -> additive discrete-time hazard):
  eqc_itt_A0.long$hazard_pos <- predict(fit2, newdata=eqc_itt_A0.long)
  eqc_itt_A1.long$hazard_pos <- predict(fit2, newdata=eqc_itt_A1.long)
  eqc_itt_A0.long$hazard_neg <- predict(fit1, newdata=eqc_itt_A0.long)
  eqc_itt_A1.long$hazard_neg <- predict(fit1, newdata=eqc_itt_A1.long)
  ### identity-link hazards are not range-restricted; clamp to [0,1]
  eqc_itt_A0.long$hazard_pos <- pmin(pmax(eqc_itt_A0.long$hazard_pos, 0), 1)
  eqc_itt_A1.long$hazard_pos <- pmin(pmax(eqc_itt_A1.long$hazard_pos, 0), 1)
  eqc_itt_A0.long$hazard_neg <- pmin(pmax(eqc_itt_A0.long$hazard_neg, 0), 1)
  eqc_itt_A1.long$hazard_neg <- pmin(pmax(eqc_itt_A1.long$hazard_neg, 0), 1)
  ### Corrected hazards under no treatment -- ADDITIVE de-biasing:
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
  eqc_itt_A0.long <- eqc_itt_A0.long[order(eqc_itt_A0.long$bootid_mrn, eqc_itt_A0.long$time_end),]
  eqc_itt_A1.long <- eqc_itt_A1.long[order(eqc_itt_A1.long$bootid_mrn, eqc_itt_A1.long$time_end),]

  ### Calculate the cumulative survival

  # lag P(no event pos)

  eqc_itt_A0.long <- eqc_itt_A0.long |>
    arrange(bootid_mrn, time_end) |>
    group_by(bootid_mrn) |>
    mutate(pnoevent_pos_lag = lag(pnoevent_pos, n=1, default=1),
           pnoevent_pos_c_lag = lag(pnoevent_pos_c, n=1, default=1)) |>
    ungroup()

  eqc_itt_A1.long <- eqc_itt_A1.long |>
    arrange(bootid_mrn, time_end) |>
    group_by(bootid_mrn) |>
    mutate(pnoevent_pos_lag = lag(pnoevent_pos, n=1, default=1)) |>
    ungroup()

  # product at each time (P(no event neg) * lag P(no event pos))

  eqc_itt_A0.long$surv_prod_lag <- eqc_itt_A0.long$pnoevent_neg * eqc_itt_A0.long$pnoevent_pos_lag
  eqc_itt_A1.long$surv_prod_lag <- eqc_itt_A1.long$pnoevent_neg * eqc_itt_A1.long$pnoevent_pos_lag
  eqc_itt_A0.long$surv_prod_c_lag <- eqc_itt_A0.long$pnoevent_neg_c * eqc_itt_A0.long$pnoevent_pos_c_lag

  # cumulative product within individual

  eqc_itt_A0.long$survival_pos <- ave(eqc_itt_A0.long$surv_prod_lag, eqc_itt_A0.long$bootid_mrn, FUN=cumprod)
  eqc_itt_A1.long$survival_pos <- ave(eqc_itt_A1.long$surv_prod_lag, eqc_itt_A1.long$bootid_mrn, FUN=cumprod)
  eqc_itt_A0.long$survival_pos_c <- ave(eqc_itt_A0.long$surv_prod_c_lag, eqc_itt_A0.long$bootid_mrn, FUN=cumprod)

  ### Calculate risk using CIF estimator

  # product at each time (haz pos * surv pos)

  eqc_itt_A0.long$risk_prod_pos <- eqc_itt_A0.long$hazard_pos * eqc_itt_A0.long$survival_pos
  eqc_itt_A1.long$risk_prod_pos <- eqc_itt_A1.long$hazard_pos * eqc_itt_A1.long$survival_pos
  eqc_itt_A0.long$risk_prod_pos_c <- eqc_itt_A0.long$hazard_pos_c * eqc_itt_A0.long$survival_pos_c

  # cumulative sum within individual

  eqc_itt_A0.long$risk_pos <- ave(eqc_itt_A0.long$risk_prod_pos, eqc_itt_A0.long$bootid_mrn, FUN=cumsum)
  eqc_itt_A1.long$risk_pos <- ave(eqc_itt_A1.long$risk_prod_pos, eqc_itt_A1.long$bootid_mrn, FUN=cumsum)
  eqc_itt_A0.long$risk_pos_c <- ave(eqc_itt_A0.long$risk_prod_pos_c, eqc_itt_A0.long$bootid_mrn, FUN=cumsum)

  # Calculate the average risk at each time point

  eqc_itt_A0.long.res <- aggregate(risk_pos ~ time_end, data=eqc_itt_A0.long, FUN=mean)
  eqc_itt_A1.long.res <- aggregate(risk_pos ~ time_end, data=eqc_itt_A1.long, FUN=mean)
  eqc_itt_A0.long.res.c <- aggregate(risk_pos_c ~ time_end, data=eqc_itt_A0.long, FUN=mean)

  message("Finished bootstrap ", i, " in ",
          round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), " minutes")

  cbind(sim = i,
        time_end = eqc_itt_A0.long.res$time_end,
        risk0 = eqc_itt_A0.long.res$risk_pos,
        risk0corr = eqc_itt_A0.long.res.c$risk_pos_c,
        risk1 = eqc_itt_A1.long.res$risk_pos)
})

boot.long <- bind_rows(lapply(boot.results, as.data.frame))
saveRDS(boot.long, paste0(res_path, "eqc.add.itt.boot.long.rds"))
