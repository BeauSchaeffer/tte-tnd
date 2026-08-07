##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Equi-Confounding Pooled ADDITIVE ITT Bootstrap
##----- Intention-to-treat
##----- last updated 2026-08-07
##-----
##----- ITT twin of equi_confounding_pooled_pp_additive_boot.R.


# Packages ----------------------------------------------------------------


library(tidyverse)
library(data.table)
library(speedglm)
library(splines)


# Data --------------------------------------------------------------------


data_Y3 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch.3/data_Y3_weekmatch.rds")
dat <- data_Y3
setDT(dat)

res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results_itt.3/"


# Boot setup --------------------------------------------------------------


num.boot <- 200

set.seed(1155)
seed <- floor(runif(num.boot)*10^8)

setkey(dat, subclass)
subclasses <- dat[, unique(subclass)]
n_sub <- length(subclasses)

ref_A1 <- data.frame(time_end = seq(1,53,1), treatment = 1,
                     sex_admin = factor("F"), age_years = 0, bmi = 0,
                     race = factor("White"), charlson_cat_fac = factor("0"),
                     ndi = 0, prior_inf = 0, tests_count = 0,
                     service_region = factor("Central valley"),
                     last_vax_infect_weeks = 0, flu_vax = 0)
ref_A0 <- ref_A1; ref_A0$treatment <- 0

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

  ## de-biased additive effect: beta_2A(t) = delta2(t) - delta1(t)
  delta1 <- predict(fit1, newdata = ref_A1) - predict(fit1, newdata = ref_A0)
  delta2 <- predict(fit2, newdata = ref_A1) - predict(fit2, newdata = ref_A0)
  time_df <- data.frame(time_end = seq(1,53,1), beta2A = delta2 - delta1)

  ## g-computation frame
  dat.boot$gmaxt <- 53
  g <- dat.boot[rep(1:nrow(dat.boot), dat.boot$gmaxt),]
  g$time_start <- ave(g$bootid_mrn, g$bootid_mrn, FUN=seq_along)
  g$time_start <- (g$time_start-1)*time_unit
  g$time_end   <- g$time_start + time_unit
  g$treatment_obs <- g$treatment
  g$lambda1 <- pmin(pmax(predict(fit1, newdata = g), 0), 1)
  g$lambda2 <- pmin(pmax(predict(fit2, newdata = g), 0), 1)
  g <- left_join(g, time_df, by = "time_end")

  g <- g |>
    arrange(bootid_mrn, time_end) |>
    group_by(bootid_mrn) |>
    mutate(
      hz0 = pmin(pmax(lambda2 + beta2A * (0 - treatment_obs), 0), 1),
      hz1 = pmin(pmax(lambda2 + beta2A * (1 - treatment_obs), 0), 1),
      surv0 = exp(-cumsum(pmin(pmax(lambda1 + hz0, 0), 1))),
      surv1 = exp(-cumsum(pmin(pmax(lambda1 + hz1, 0), 1))),
      risk0 = cumsum(hz0 * lag(surv0, default = 1)),
      risk1 = cumsum(hz1 * lag(surv1, default = 1))
    ) |>
    ungroup()

  res <- g |> group_by(time_end) |> summarise(risk0 = mean(risk0), risk1 = mean(risk1), .groups = "drop")

  message("Finished bootstrap ", i, " in ",
          round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), " minutes")

  cbind(sim = i, time_end = res$time_end, risk0 = res$risk0, risk1 = res$risk1)
})

boot.long <- bind_rows(lapply(boot.results, as.data.frame))
saveRDS(boot.long, paste0(res_path, "eqc.add.itt.boot.long.rds"))
