##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Equi-Confounding Analysis Pooled
##----- Per-protocol, no censoring weights
##----- last updated 2026-07-28


# Packages ----------------------------------------------------------------


library(tidyverse)
library(data.table)
library(speedglm)
library(splines)


# Data --------------------------------------------------------------------


data_Y3 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch.3/data_Y3_weekmatch.rds")
dat <- data_Y3
setDT(dat)

res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results_pp.3/"


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


eqc_pooled_pp_fit1 <- speedglm(Y_neg ~ ns(time_end, knots = c(10,20,30,40,50))*treatment +
                                  # demographic
                                  sex_admin + age_years + bmi + race + charlson_cat_fac +
                                  # other
                                  ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks + 
                                  # NEC
                                  flu_vax,
                                data=dat.long.pp,
                                family=binomial())

saveRDS(eqc_pooled_pp_fit1, paste0(res_path,"eqc_pooled_pp_fit1.rds"))



eqc_pooled_pp_fit2 <- speedglm(Y_pos ~ ns(time_end, knots = c(10,20,30,40,50))*treatment +
                                  # demographic
                                  sex_admin + age_years + bmi + race + charlson_cat_fac +
                                  # other
                                  ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks + 
                                  # NEC
                                  flu_vax,
                                data=dat.long.pp,
                                family=binomial())

saveRDS(eqc_pooled_pp_fit2, paste0(res_path,"eqc_pooled_pp_fit2.rds")) 


# PP Survival and Risk ---------------------------------------------------


dat$gmaxt <- 53

### G formula data setup A=0
eqc_pp_A0.long <- dat[rep(1:nrow(dat), dat$gmaxt),]
eqc_pp_A0.long$time_start <- ave(eqc_pp_A0.long$fake_mrn, eqc_pp_A0.long$fake_mrn, FUN=seq_along)
eqc_pp_A0.long$time_start <- (eqc_pp_A0.long$time_start-1)*time_unit
eqc_pp_A0.long$time_end <- eqc_pp_A0.long$time_start+time_unit
eqc_pp_A0.long$treatment <- 0

### G formula data setup A=1
eqc_pp_A1.long <- dat[rep(1:nrow(dat), dat$gmaxt),]
eqc_pp_A1.long$time_start <- ave(eqc_pp_A1.long$fake_mrn, eqc_pp_A1.long$fake_mrn, FUN=seq_along)
eqc_pp_A1.long$time_start <- (eqc_pp_A1.long$time_start-1)*time_unit
eqc_pp_A1.long$time_end <- eqc_pp_A1.long$time_start+time_unit
eqc_pp_A1.long$treatment <- 1

### Calculate predicted hazards:
eqc_pp_A0.long$hazard_pos <- predict(eqc_pooled_pp_fit2, newdata=eqc_pp_A0.long, type="response")
eqc_pp_A1.long$hazard_pos <- predict(eqc_pooled_pp_fit2, newdata=eqc_pp_A1.long, type="response")
eqc_pp_A0.long$hazard_neg <- predict(eqc_pooled_pp_fit1, newdata=eqc_pp_A0.long, type="response")
eqc_pp_A1.long$hazard_neg <- predict(eqc_pooled_pp_fit1, newdata=eqc_pp_A1.long, type="response")
### Corrected hazards under no treatment
eqc_pp_A0.long$hazard_pos_c <- eqc_pp_A0.long$hazard_pos * (eqc_pp_A1.long$hazard_neg / eqc_pp_A0.long$hazard_neg)

### Calculate (1 - hazard)
eqc_pp_A0.long$pnoevent_pos <- 1 - eqc_pp_A0.long$hazard_pos
eqc_pp_A1.long$pnoevent_pos <- 1 - eqc_pp_A1.long$hazard_pos
eqc_pp_A0.long$pnoevent_neg <- 1 - eqc_pp_A0.long$hazard_neg
eqc_pp_A1.long$pnoevent_neg <- 1 - eqc_pp_A1.long$hazard_neg
### Corrected (1 - hazard) under no treatment
eqc_pp_A0.long$pnoevent_pos_c <- 1 - eqc_pp_A0.long$hazard_pos_c

### Sort the data by ID, time
eqc_pp_A0.long <- eqc_pp_A0.long[order(eqc_pp_A0.long$fake_mrn, eqc_pp_A0.long$time_end),] 
eqc_pp_A1.long <- eqc_pp_A1.long[order(eqc_pp_A1.long$fake_mrn, eqc_pp_A1.long$time_end),] 

### Calculate the cumulative survival 

# lag P(no event pos)

eqc_pp_A0.long <- eqc_pp_A0.long |> 
  arrange(fake_mrn, time_end) |> 
  group_by(fake_mrn) |> 
  mutate(pnoevent_pos_lag = lag(pnoevent_pos, n=1, default=1),
         pnoevent_pos_c_lag = lag(pnoevent_pos_c, n=1, default=1)) |> 
  ungroup()

eqc_pp_A1.long <- eqc_pp_A1.long |> 
  arrange(fake_mrn, time_end) |> 
  group_by(fake_mrn) |> 
  mutate(pnoevent_pos_lag = lag(pnoevent_pos, n=1, default=1)) |> 
  ungroup()

# product at each time (P(no event neg) * lag P(no event pos))

eqc_pp_A0.long$surv_prod_lag <- eqc_pp_A0.long$pnoevent_neg * eqc_pp_A0.long$pnoevent_pos_lag
eqc_pp_A1.long$surv_prod_lag <- eqc_pp_A1.long$pnoevent_neg * eqc_pp_A1.long$pnoevent_pos_lag
eqc_pp_A0.long$surv_prod_c_lag <- eqc_pp_A0.long$pnoevent_neg * eqc_pp_A0.long$pnoevent_pos_c_lag

# cumulative product within individual

eqc_pp_A0.long$survival_pos <- ave(eqc_pp_A0.long$surv_prod_lag, eqc_pp_A0.long$fake_mrn, FUN=cumprod)
eqc_pp_A1.long$survival_pos <- ave(eqc_pp_A1.long$surv_prod_lag, eqc_pp_A1.long$fake_mrn, FUN=cumprod)
eqc_pp_A0.long$survival_pos_c <- ave(eqc_pp_A0.long$surv_prod_c_lag, eqc_pp_A0.long$fake_mrn, FUN=cumprod)

### Calculate risk using CIF estimator

# product at each time (haz pos * surv pos)

eqc_pp_A0.long$risk_prod_pos <- eqc_pp_A0.long$hazard_pos * eqc_pp_A0.long$survival_pos
eqc_pp_A1.long$risk_prod_pos <- eqc_pp_A1.long$hazard_pos * eqc_pp_A1.long$survival_pos
eqc_pp_A0.long$risk_prod_pos_c <- eqc_pp_A0.long$hazard_pos_c * eqc_pp_A0.long$survival_pos_c

# cumulative sum within individual

eqc_pp_A0.long$risk_pos <- ave(eqc_pp_A0.long$risk_prod_pos, eqc_pp_A0.long$fake_mrn, FUN=cumsum)
eqc_pp_A1.long$risk_pos <- ave(eqc_pp_A1.long$risk_prod_pos, eqc_pp_A1.long$fake_mrn, FUN=cumsum)
eqc_pp_A0.long$risk_pos_c <- ave(eqc_pp_A0.long$risk_prod_pos_c, eqc_pp_A0.long$fake_mrn, FUN=cumsum)

# Calculate the average risk at each time point

eqc_pp_A0.long.res <- aggregate(risk_pos ~ time_end, data=eqc_pp_A0.long, FUN=mean)
eqc_pp_A1.long.res <- aggregate(risk_pos ~ time_end, data=eqc_pp_A1.long, FUN=mean)
eqc_pp_A0.long.res.c <- aggregate(risk_pos_c ~ time_end, data=eqc_pp_A0.long, FUN=mean)


# save point estimate risk curves in bootstrap-compatible format
eqc.pp.risk.pointest <- tibble(
  sim = 0L,  # 0 = main analysis (bootstraps are 1..B)
  time_end = eqc_pp_A0.long.res$time_end,
  risk0 = eqc_pp_A0.long.res$risk_pos,
  risk0corr = eqc_pp_A0.long.res.c$risk_pos_c,
  risk1 = eqc_pp_A1.long.res$risk_pos
)

saveRDS(eqc.pp.risk.pointest, paste0(res_path, "eqc.pp.risk.pointest.rds"))

