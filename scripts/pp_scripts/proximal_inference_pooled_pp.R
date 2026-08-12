##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Proximal Inference Analysis Pooled
##----- Per-protocol, no censoring weights
##----- last updated 2026-08-07


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


# Downsample --------------------------------------------------------------


# subclass_ids <- data_Y3 |> dplyr::select(subclass) |> unique()
# set.seed(345)
# subclass_ids_subset <- dplyr::slice_sample(subclass_ids, n=10000)
# dat_downsamp <- data_Y3 |> dplyr::filter(subclass %in% subclass_ids_subset$subclass) |> droplevels()
# rm(subclass_ids, subclass_ids_subset)


# PP long format expansion -----------------------------------------------


### calc number of rows needed for each individual
time_unit <- 1

### ensure at least 1 row for each individual
dat$max_units <- ceiling(dat$Y3_pp_t_trunc/time_unit)+1
dat.long.pp <- dat[rep(1:nrow(dat), dat$max_units),]

### variable that represents the start and end time corresponding to each row of observation
dat.long.pp$time_start <- ave(dat.long.pp$fake_mrn, dat.long.pp$fake_mrn, FUN=seq_along)
dat.long.pp$time_start <- (dat.long.pp$time_start-1)*time_unit
dat.long.pp$time_end <- dat.long.pp$time_start+time_unit

### modify the Y and C variables so that they are only equal to 1 if the 
### event/censoring happened in that time interval
dat.long.pp$Y_pos <- ifelse(
  dat.long.pp$Y3_pp_trunc == 2 &
    dat.long.pp$Y3_pp_t_trunc == dat.long.pp$time_start,
  1, 0
)

dat.long.pp$Y_neg <- ifelse(
  dat.long.pp$Y3_pp_trunc == 1 &
    dat.long.pp$Y3_pp_t_trunc == dat.long.pp$time_start,
  1, 0
)

dat.long.pp$C <- ifelse(
  dat.long.pp$Y3_pp_trunc == 0 &
    dat.long.pp$Y3_pp_t_trunc == dat.long.pp$time_start,
  1, 0
)


dat.long.pp$Y_pos <- ifelse(dat.long.pp$C==1, NA, dat.long.pp$Y_pos)
dat.long.pp$Y_neg <- ifelse(dat.long.pp$C==1, NA, dat.long.pp$Y_neg)


# PP Pooled Logistic -----------------------------------------------------

### time interacting with all variables, note ns()*()
### mem pressure peaks around 75 GB
### 100 GB with speedglm - works in a few mins

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
saveRDS(prox_pooled_pp_s1, paste0(res_path,"prox_pooled_pp_s1.rds"))

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
saveRDS(prox_pooled_pp_s2, paste0(res_path,"prox_pooled_pp_s2.rds"))

# # sanity check against cox model
# # stage 2 with no tx interaction
# prox_pooled_pp_s2_noint <- speedglm(Y_pos ~ ns(time_end, knots = c(10,20,30,40,50)) + treatment +
#                                 # demographic
#                                 sex_admin + age_years + bmi + race + charlson_cat_fac +
#                                 # other
#                                 ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
#                                 # predictions from stage 1
#                                 p_pp,
#                               # no NEC
#                               data=dat.long.pp,
#                               family=binomial(),
#                               sparse = FALSE)
# saveRDS(prox_pooled_pp_s2_noint, paste0(res_path,"prox_pooled_pp_s2_noint.rds"))
# exp(prox_pooled_pp_s2_noint$coefficients["treatment"])

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
saveRDS(prox_pooled_pp_obs, paste0(res_path,"prox_pooled_pp_obs.rds"))


# PP Survival and Risk ---------------------------------------------------

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

### Evaluate BOTH anchor hazards at the observed treatment: treatment has no causal
### effect on the test-negative event (NCO), and the test-positive counterfactual is
### applied via the switching function below; anchoring at the intervention value
### would double-count the treatment effect in that step.
prox_pp_A0.long$treatment <- prox_pp_A0.long$treatment_obs
prox_pp_A1.long$treatment <- prox_pp_A1.long$treatment_obs
prox_pp_A0.long$hazard_neg <- predict(prox_pooled_pp_s1, newdata=prox_pp_A0.long, type="response") 
prox_pp_A1.long$hazard_neg <- predict(prox_pooled_pp_s1, newdata=prox_pp_A1.long, type="response")
# ### predicted hazards testing POSITIVE from stage 2 model
### overwrite below
# prox_pp_A0.long$hazard_pos <- predict(prox_pooled_pp_s2, newdata=prox_pp_A0.long, type="response")
# prox_pp_A1.long$hazard_pos <- predict(prox_pooled_pp_s2, newdata=prox_pp_A1.long, type="response")

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

### calculate (1 - hazard NEGATIVE)
prox_pp_A0.long$pnoevent_neg <- 1 - prox_pp_A0.long$hazard_neg
prox_pp_A1.long$pnoevent_neg <- 1 - prox_pp_A1.long$hazard_neg

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

### lag (1 - hazard NEGATIVE)
prox_pp_A0.long <- prox_pp_A0.long |> 
  arrange(fake_mrn, time_end) |> 
  group_by(fake_mrn) |> 
  mutate(pnoevent_neg_lag = lag(pnoevent_neg, n=1, default=1)) |> 
  ungroup()

prox_pp_A1.long <- prox_pp_A1.long |> 
  arrange(fake_mrn, time_end) |> 
  group_by(fake_mrn) |> 
  mutate(pnoevent_neg_lag = lag(pnoevent_neg, n=1, default=1)) |> 
  ungroup()


### Aalen–Johansen estimator

prox_pp_A0.long$surv_prod <- prox_pp_A0.long$pnoevent_neg * prox_pp_A0.long$pnoevent_pos_lag 
prox_pp_A1.long$surv_prod <- prox_pp_A1.long$pnoevent_neg * prox_pp_A1.long$pnoevent_pos_lag

prox_pp_A0.long$survival <- ave(prox_pp_A0.long$surv_prod, prox_pp_A0.long$fake_mrn, FUN=cumprod)
prox_pp_A1.long$survival <- ave(prox_pp_A1.long$surv_prod, prox_pp_A1.long$fake_mrn, FUN=cumprod)

prox_pp_A0.long$risk_prod <- prox_pp_A0.long$hazard_pos * prox_pp_A0.long$survival
prox_pp_A1.long$risk_prod <- prox_pp_A1.long$hazard_pos * prox_pp_A1.long$survival

prox_pp_A0.long$risk_pos <- ave(prox_pp_A0.long$risk_prod, prox_pp_A0.long$fake_mrn, FUN=cumsum)
prox_pp_A1.long$risk_pos <- ave(prox_pp_A1.long$risk_prod, prox_pp_A1.long$fake_mrn, FUN=cumsum)

prox_pp_A0.long.res <- aggregate(risk_pos ~ time_end, data=prox_pp_A0.long, FUN=mean)
prox_pp_A1.long.res <- aggregate(risk_pos ~ time_end, data=prox_pp_A1.long, FUN=mean)

### save point estimate risk curves in bootstrap-compatible format
pci.pp.risk.pointest <- tibble(
  sim = 0L,  # 0 = main analysis (bootstraps are 1..B)
  time_end = prox_pp_A0.long.res$time_end,
  risk0 = prox_pp_A0.long.res$risk_pos,
  risk1 = prox_pp_A1.long.res$risk_pos
)

saveRDS(pci.pp.risk.pointest, paste0(res_path, "pci.pp.risk.pointest.rds")) 

