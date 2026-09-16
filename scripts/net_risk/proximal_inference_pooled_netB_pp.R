##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Proximal Inference Analysis Pooled -- NET RISK, Option B
##----- Per-protocol, no censoring weights
##----- last updated 2026-09-14


# Packages ----------------------------------------------------------------


library(tidyverse)
library(data.table)
library(speedglm)
library(splines)


# Data --------------------------------------------------------------------


### Net risk replaces the single three-level Y3 dataset with two outcomes
### defined on ONE risk set, {T2 > t}:
###   Y_pos = terminal indicator, 1 in the week of the first positive
###   Y_neg = recurrent indicator, 1 in every week a negative test occurred
### Under Option B both models are fit on the same person-week rows.

data_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch.4/"

data_Y2  <- read_rds(paste0(data_path, "data_Y2_weekmatch.rds"))
neg_hist <- read_rds(paste0(data_path, "neg_hist_weekmatch.rds"))

dat <- data_Y2
setDT(dat)

res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results_pp.5/"
dir.create(res_path, recursive = TRUE, showWarnings = FALSE)
stopifnot(dir.exists(res_path))

### Input checks
mod_vars <- c("treatment",
              "sex_admin", "age_years", "bmi", "race", "charlson_cat_fac",
              "ndi", "prior_inf", "tests_count", "service_region",
              "last_vax_infect_weeks",
              "flu_vax")

stopifnot(!anyDuplicated(dat$fake_mrn))
stopifnot(!anyNA(dat |> dplyr::select(all_of(mod_vars))))
stopifnot(!anyNA(dat$Y2_pp_t_trunc), !anyNA(dat$Y2_pp_trunc))
stopifnot(all(dat$Y2_pp_t_trunc > 0))
stopifnot(all(neg_hist$fake_mrn %in% dat$fake_mrn))

message("Cohort: ", nrow(dat), " individuals, ",
        n_distinct(dat$subclass), " matched pairs")
message("Positive events (PP, truncated): ", sum(dat$Y2_pp_trunc))


# Downsample --------------------------------------------------------------


# subclass_ids <- data_Y2 |> dplyr::select(subclass) |> unique()
# set.seed(345)
# subclass_ids_subset <- dplyr::slice_sample(subclass_ids, n=10000)
# dat_downsamp <- data_Y2 |> dplyr::filter(subclass %in% subclass_ids_subset$subclass) |> droplevels()
# rm(subclass_ids, subclass_ids_subset)


# PP long format expansion -----------------------------------------------


### calc number of rows needed for each individual
time_unit <- 1

### ensure at least 1 row for each individual
dat$max_units <- ceiling(dat$Y2_pp_t_trunc/time_unit)+1
dat.long.pp <- dat[rep(1:nrow(dat), dat$max_units),]

### variable that represents the start and end time corresponding to each row of observation
dat.long.pp$time_start <- ave(dat.long.pp$fake_mrn, dat.long.pp$fake_mrn, FUN=seq_along)
dat.long.pp$time_start <- (dat.long.pp$time_start-1)*time_unit
dat.long.pp$time_end <- dat.long.pp$time_start+time_unit

### modify the Y and C variables so that they are only equal to 1 if the 
### event/censoring happened in that time interval
dat.long.pp$Y_pos <- ifelse(
  dat.long.pp$Y2_pp_trunc == 1 &
    dat.long.pp$Y2_pp_t_trunc == dat.long.pp$time_start,
  1, 0
)

dat.long.pp$C <- ifelse(
  dat.long.pp$Y2_pp_trunc == 0 &
    dat.long.pp$Y2_pp_t_trunc == dat.long.pp$time_start,
  1, 0
)

### Y_neg is recurrent: 1 in every week a negative test occurred, marked on
### the same rows as Y_pos. Negatives at the cap are kept (<=) since time is
### in whole weeks; same-week repeats collapse to a single indicator.

neg_wk <- neg_hist |>
  inner_join(dat |> dplyr::select(fake_mrn, cap = Y2_pp_t_trunc), by = "fake_mrn") |>
  filter(neg_t <= cap) |>
  distinct(fake_mrn, time_start = neg_t)

setDT(neg_wk)
setDT(dat.long.pp)

dat.long.pp[, Y_neg := 0]
dat.long.pp[neg_wk, on = c("fake_mrn", "time_start"), Y_neg := 1]

dat.long.pp$Y_pos <- ifelse(dat.long.pp$C==1, NA, dat.long.pp$Y_pos)
dat.long.pp$Y_neg <- ifelse(dat.long.pp$C==1, NA, dat.long.pp$Y_neg)

message("Person-weeks: ", nrow(dat.long.pp),
        " (", sum(dat.long.pp$C==1), " dropped as censoring rows)")
message("Positive person-weeks: ", sum(dat.long.pp$Y_pos, na.rm = TRUE),
        " | negative person-weeks: ", sum(dat.long.pp$Y_neg, na.rm = TRUE))


# PP Pooled Logistic -----------------------------------------------------

### time interacting with all variables, note ns()*()
### mem pressure peaks around 75 GB
### 100 GB with speedglm - works in a few mins

### stage 1 is now the recurrent negative-testing hazard on {T2 > t} rather
### than a cause-specific first-test hazard. Specification unchanged.

prox_pooled_pp_s1 <- speedglm(Y_neg ~ ns(time_end, knots = c(10,20,30,40,50))*(treatment +
                                                                            # demographic
                                                                            sex_admin + age_years + bmi + race + charlson_cat_fac +
                                                                            # other
                                                                            ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks + 
                                                                            # NEC
                                                                            flu_vax),
                               data=dat.long.pp,
                               family=binomial(),
                              sparse = FALSE)
saveRDS(prox_pooled_pp_s1, paste0(res_path,"netB.prox_pooled_pp_s1.rds"))

dat.long.pp$p_pp <- predict(prox_pooled_pp_s1, newdata = dat.long.pp)

prox_pooled_pp_s2 <- speedglm(Y_pos ~ ns(time_end, knots = c(10,20,30,40,50))*treatment +
                                 # demographic
                                 sex_admin + age_years + bmi + race + charlson_cat_fac +
                                 # other
                                 ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks + 
                                 # predictions from stage 1
                                 p_pp,
                               # no NEC
                               data=dat.long.pp,
                               family=binomial(),
                              sparse = FALSE)
saveRDS(prox_pooled_pp_s2, paste0(res_path,"netB.prox_pooled_pp_s2.rds"))

prox_pooled_pp_obs <- speedglm(Y_pos ~ ns(time_end, knots = c(10,20,30,40,50))*treatment +
                                  # demographic
                                  sex_admin + age_years + bmi + race + charlson_cat_fac +
                                  # other
                                  ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks + 
                                  # NEC
                                  flu_vax,
                                data=dat.long.pp,
                                family=binomial(),
                               sparse = FALSE)
saveRDS(prox_pooled_pp_obs, paste0(res_path,"netB.prox_pooled_pp_obs.rds"))


# PP Survival and Risk ---------------------------------------------------

### Net risk is a SINGLE decrement: only the positive hazard removes anyone
### from the risk set, so survival is prod(1 - hazard_pos) alone. The negative
### hazard does not enter the survival product -- under PCI it enters only
### through the stage 1 bridge p_pp.

dat$gmaxt <- 53

### G formula data setup A=0
prox_pp_A0.long <- dat[rep(1:nrow(dat), dat$gmaxt),]
prox_pp_A0.long$time_start <- ave(prox_pp_A0.long$fake_mrn, prox_pp_A0.long$fake_mrn, FUN=seq_along)
prox_pp_A0.long$time_start <- (prox_pp_A0.long$time_start-1)*time_unit
prox_pp_A0.long$time_end <- prox_pp_A0.long$time_start+time_unit
prox_pp_A0.long$treatment_obs <- prox_pp_A0.long$treatment
prox_pp_A0.long$treatment <- 0

### G formula data setup A=1
prox_pp_A1.long <- dat[rep(1:nrow(dat), dat$gmaxt),]
prox_pp_A1.long$time_start <- ave(prox_pp_A1.long$fake_mrn, prox_pp_A1.long$fake_mrn, FUN=seq_along)
prox_pp_A1.long$time_start <- (prox_pp_A1.long$time_start-1)*time_unit
prox_pp_A1.long$time_end <- prox_pp_A1.long$time_start+time_unit
prox_pp_A1.long$treatment_obs <- prox_pp_A1.long$treatment
prox_pp_A1.long$treatment <- 1

### stage 1 linear predictor under each intervention
prox_pp_A0.long$p_pp <- predict(prox_pooled_pp_s1, newdata=prox_pp_A0.long, type="link") 
prox_pp_A1.long$p_pp <- predict(prox_pooled_pp_s1, newdata=prox_pp_A1.long, type="link") 

### Evaluate the anchor hazard at the observed treatment: the test-positive
### counterfactual is applied via the switching function below, so anchoring at
### the intervention value would double-count the treatment effect in that step.
prox_pp_A0.long$treatment <- prox_pp_A0.long$treatment_obs
prox_pp_A1.long$treatment <- prox_pp_A1.long$treatment_obs

### predicted hazards testing POSITIVE from observed model (at observed treatment)
prox_pp_A0.long$hazard_pos_obs <- predict(prox_pooled_pp_obs, newdata=prox_pp_A0.long, type="response")
prox_pp_A1.long$hazard_pos_obs <- predict(prox_pooled_pp_obs, newdata=prox_pp_A1.long, type="response")

### referent data frames for extracting stage 2 treatment contrasts at each time t

df_ref_A1 <- data.frame(time_end=seq(1,53,1),
                        treatment=1,
                        sex_admin=factor("F"),
                        age_years=0,
                        bmi=0,
                        race=factor("White"),
                        charlson_cat_fac=factor("0"),
                        ndi=0,
                        prior_inf=0,
                        tests_count=0,
                        service_region=factor("Central valley"),
                        last_vax_infect_weeks=0,
                        p_pp=0)

df_ref_A0 <- data.frame(time_end=seq(1,53,1),
                        treatment=0,
                        sex_admin=factor("F"),
                        age_years=0,
                        bmi=0,
                        race=factor("White"),
                        charlson_cat_fac=factor("0"),
                        ndi=0,
                        prior_inf=0,
                        tests_count=0,
                        service_region=factor("Central valley"),
                        last_vax_infect_weeks=0,
                        p_pp=0)

haz_ref_A1 <- predict(prox_pooled_pp_s2, newdata=df_ref_A1, type = "link")
haz_ref_A0 <- predict(prox_pooled_pp_s2, newdata=df_ref_A0, type = "link")

time_df <- data.frame(time_end=seq(1,53,1),
                      logHR=haz_ref_A1-haz_ref_A0)

prox_pp_A0.long <- left_join(prox_pp_A0.long, time_df, by="time_end")
prox_pp_A1.long <- left_join(prox_pp_A1.long, time_df, by="time_end")

### switching function

### removing treatment from treated
### take hazard from untreated, remove treatment from treated
### negative log HR to remove
prox_pp_A0.long$hazard_pos <- prox_pp_A0.long$hazard_pos_obs * exp(-prox_pp_A0.long$logHR * prox_pp_A0.long$treatment_obs)
### adding treated to untreated
### take hazard from treated, add treatment to untreated
### positive log HR to add
prox_pp_A1.long$hazard_pos <- prox_pp_A1.long$hazard_pos_obs * exp(prox_pp_A1.long$logHR * (1-prox_pp_A1.long$treatment_obs))

### compute survival and cumulative incidence

### calculate (1 - hazard POSITIVE)
prox_pp_A0.long$pnoevent_pos <- 1 - prox_pp_A0.long$hazard_pos
prox_pp_A1.long$pnoevent_pos <- 1 - prox_pp_A1.long$hazard_pos

### sort the data by ID, time
prox_pp_A0.long <- prox_pp_A0.long[order(prox_pp_A0.long$fake_mrn, prox_pp_A0.long$time_end),] 
prox_pp_A1.long <- prox_pp_A1.long[order(prox_pp_A1.long$fake_mrn, prox_pp_A1.long$time_end),]

### lag (1 - hazard POSITIVE)
prox_pp_A0.long <- prox_pp_A0.long |> 
  arrange(fake_mrn, time_end) |> 
  group_by(fake_mrn) |> 
  mutate(pnoevent_pos_lag = lag(pnoevent_pos, n=1, default=1)) |> 
  ungroup()

prox_pp_A1.long <- prox_pp_A1.long |> 
  arrange(fake_mrn, time_end) |> 
  group_by(fake_mrn) |> 
  mutate(pnoevent_pos_lag = lag(pnoevent_pos, n=1, default=1)) |> 
  ungroup()


### Kaplan-Meier type estimator (single decrement -- lagged positive hazard
### only, no competing test-negative term)

prox_pp_A0.long$survival <- ave(prox_pp_A0.long$pnoevent_pos_lag, prox_pp_A0.long$fake_mrn, FUN=cumprod)
prox_pp_A1.long$survival <- ave(prox_pp_A1.long$pnoevent_pos_lag, prox_pp_A1.long$fake_mrn, FUN=cumprod)

prox_pp_A0.long$risk_prod <- prox_pp_A0.long$hazard_pos * prox_pp_A0.long$survival
prox_pp_A1.long$risk_prod <- prox_pp_A1.long$hazard_pos * prox_pp_A1.long$survival

prox_pp_A0.long$risk_pos <- ave(prox_pp_A0.long$risk_prod, prox_pp_A0.long$fake_mrn, FUN=cumsum)
prox_pp_A1.long$risk_pos <- ave(prox_pp_A1.long$risk_prod, prox_pp_A1.long$fake_mrn, FUN=cumsum)

prox_pp_A0.long.res <- aggregate(risk_pos ~ time_end, data=prox_pp_A0.long, FUN=mean)
prox_pp_A1.long.res <- aggregate(risk_pos ~ time_end, data=prox_pp_A1.long, FUN=mean)

### save point estimate risk curves in bootstrap-compatible format
netB.pci.pp.risk.pointest <- tibble(
  sim = 0L,  # 0 = main analysis (bootstraps are 1..B)
  time_end = prox_pp_A0.long.res$time_end,
  risk0 = prox_pp_A0.long.res$risk_pos,
  risk1 = prox_pp_A1.long.res$risk_pos
)

saveRDS(netB.pci.pp.risk.pointest, paste0(res_path, "netB.pci.pp.risk.pointest.rds")) 
