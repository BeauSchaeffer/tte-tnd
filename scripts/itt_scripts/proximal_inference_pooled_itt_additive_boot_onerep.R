##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Proximal Inference Pooled ADDITIVE ITT Bootstrap
##----- Intention-to-treat
##----- ** SINGLE BOOTSTRAP REPLICATE FOR USE WITH ARRAY **
##----- last updated 2026-08-07
##-----
##----- ITT twin of proximal_inference_pooled_pp_additive_boot_onerep.R.


# Packages ----------------------------------------------------------------


library(tidyverse)
library(data.table)
library(speedglm)
library(splines)


# Data --------------------------------------------------------------------


data_Y3 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch.3/data_Y3_weekmatch.rds")
dat <- data_Y3
setDT(dat)

res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results_itt.5/pci_add_boot_reps"
dir.create(res_path, showWarnings = FALSE, recursive = TRUE)


# Boot --------------------------------------------------------------------


args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])

num.boot <- 200
set.seed(1155)
seed <- floor(runif(num.boot)*10^8)
set.seed(seed[i])

setkey(dat, subclass)
subclasses <- dat[, unique(subclass)]
n_sub <- length(subclasses)

time_unit <- 1

t0 <- Sys.time()
message("Starting PCI ITT additive bootstrap ", i, " (seed=", seed[i], ")")

## resample matched sets, unique ids per draw
samp_sub <- sample(subclasses, size = n_sub, replace = TRUE)
map <- data.table(j = seq_along(samp_sub), subclass = samp_sub)
dat.boot <- dat[map, on = "subclass", allow.cartesian = TRUE]
dat.boot[, bootid := j]
dat.boot[, bootid_mrn := .GRP, by = .(bootid, fake_mrn)]

## long format
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

## stage 1 (additive NCO model, incl. NCE)
s1 <- speedglm(Y_neg ~ ns(time_end, knots = c(10,20,30,40,50))*(treatment +
                 sex_admin + age_years + bmi + race + charlson_cat_fac +
                 ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
                 flu_vax),
               data = dat.long, family = gaussian(), sparse=FALSE)

## bridge = fitted stage-1 hazard minus time baseline (avoids stage-2 collinearity)
ref_week <- data.frame(time_end = seq(1,53,1), treatment = 0,
                       sex_admin = factor("F"), age_years = 0, bmi = 0,
                       race = factor("White"), charlson_cat_fac = factor("0"),
                       ndi = 0, prior_inf = 0, tests_count = 0,
                       service_region = factor("Central valley"),
                       last_vax_infect_weeks = 0, flu_vax = 0)
ref_week$baseline1 <- predict(s1, newdata = ref_week)
dat.long$fitted1 <- predict(s1, newdata = dat.long)
dat.long <- left_join(dat.long, ref_week[, c("time_end","baseline1")], by = "time_end")
dat.long$mu <- dat.long$fitted1 - dat.long$baseline1

## stage 2 (bridge, no NCE) and observed primary model (NCE, no bridge)
s2  <- speedglm(Y_pos ~ ns(time_end, knots = c(10,20,30,40,50))*treatment +
                  sex_admin + age_years + bmi + race + charlson_cat_fac +
                  ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
                  mu,
                data = dat.long, family = gaussian(), sparse=FALSE)
obs <- speedglm(Y_pos ~ ns(time_end, knots = c(10,20,30,40,50))*treatment +
                  sex_admin + age_years + bmi + race + charlson_cat_fac +
                  ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
                  flu_vax,
                data = dat.long, family = gaussian(), sparse=FALSE)

## de-biased additive effect beta_2A(t) from stage-2 treatment contrast (mu = 0)
df_ref_A1 <- data.frame(time_end = seq(1,53,1), treatment = 1,
                        sex_admin = factor("F"), age_years = 0, bmi = 0,
                        race = factor("White"), charlson_cat_fac = factor("0"),
                        ndi = 0, prior_inf = 0, tests_count = 0,
                        service_region = factor("Central valley"),
                        last_vax_infect_weeks = 0, mu = 0)
df_ref_A0 <- df_ref_A1; df_ref_A0$treatment <- 0
beta2A <- predict(s2, newdata = df_ref_A1) - predict(s2, newdata = df_ref_A0)
time_df <- data.frame(time_end = seq(1,53,1), beta2A = beta2A)

## g-computation frame
dat.boot$gmaxt <- 53
g <- dat.boot[rep(1:nrow(dat.boot), dat.boot$gmaxt),]
g$time_start <- ave(g$bootid_mrn, g$bootid_mrn, FUN=seq_along)
g$time_start <- (g$time_start-1)*time_unit
g$time_end   <- g$time_start + time_unit
g$treatment_obs <- g$treatment
g$lambda1 <- pmin(pmax(predict(s1,  newdata = g), 0), 1)
g$lambda2 <- pmin(pmax(predict(obs, newdata = g), 0), 1)
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

out <- cbind(sim = i, time_end = res$time_end, risk0 = res$risk0, risk1 = res$risk1)
saveRDS(out, sprintf("%s/pci_add_itt_boot_rep_%03d.rds", res_path, i))
