##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Proximal Inference Analysis Pooled PP Bootstrap
##----- Per-protocol, no censoring weights
##----- ** SINGLE BOOTSTRAP REPLICATE FOR USE WITH ARRAY **
##----- last updated 2026-08-04


# Packages ----------------------------------------------------------------


library(tidyverse)
library(data.table)
library(speedglm)
library(splines)


# Data --------------------------------------------------------------------


data_Y3 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch.3/data_Y3_weekmatch.rds")
dat <- data_Y3
setDT(dat)

res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results_pp.4/pci_boot_reps"
dir.create(res_path, showWarnings = FALSE, recursive = TRUE)


# Boot --------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])

num.boot <- 200

set.seed(1155)
seed <- floor(runif(num.boot)*10^8)
set.seed(seed[i])

setDT(dat)
setkey(dat, subclass)
subclasses <- dat[, unique(subclass)]
n_sub <- length(subclasses)

t0 <- Sys.time()

message(
  "Starting PCI PP pooled bootstrap ", i,
  " (seed=", seed[i], ") at ",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)

# select matched pairs
samp_sub <- sample(subclasses, size = n_sub, replace = TRUE)

# build boot dataset efficiently via one join:
# map draw index j -> sampled subclass, then join to replicate all rows per subclass
map <- data.table(j = seq_along(samp_sub), subclass = samp_sub)
dat.boot <- dat[map, on = "subclass", allow.cartesian = TRUE]
# new cluster/matched pair id per draw
dat.boot[, bootid := j]
# new individual id per draw
# use data.table special group index variable
dat.boot[, bootid_mrn := .GRP, by = .(bootid, fake_mrn)]

# long format data
time_unit <- 1

dat.boot$max_units <- ceiling(dat.boot$Y3_pp_t_trunc/time_unit)+1
dat.long.boot.pp <- dat.boot[rep(1:nrow(dat.boot), dat.boot$max_units),]

dat.long.boot.pp$time_start <- ave(dat.long.boot.pp$bootid_mrn, dat.long.boot.pp$bootid_mrn, FUN=seq_along)
dat.long.boot.pp$time_start <- (dat.long.boot.pp$time_start-1)*time_unit
dat.long.boot.pp$time_end <- dat.long.boot.pp$time_start+time_unit

# recommended add
dat.long.boot.pp <- dat.long.boot.pp[order(dat.long.boot.pp$bootid_mrn, dat.long.boot.pp$time_end),]

dat.long.boot.pp$Y_pos <- ifelse(
  dat.long.boot.pp$Y3_pp_trunc == 2 &
    dat.long.boot.pp$Y3_pp_t_trunc == dat.long.boot.pp$time_start,
  1, 0
)

dat.long.boot.pp$Y_neg <- ifelse(
  dat.long.boot.pp$Y3_pp_trunc == 1 &
    dat.long.boot.pp$Y3_pp_t_trunc == dat.long.boot.pp$time_start,
  1, 0
)

dat.long.boot.pp$C <- ifelse(
  dat.long.boot.pp$Y3_pp_trunc == 0 &
    dat.long.boot.pp$Y3_pp_t_trunc == dat.long.boot.pp$time_start,
  1, 0
)

dat.long.boot.pp$Y_pos <- ifelse(dat.long.boot.pp$C==1, NA, dat.long.boot.pp$Y_pos)
dat.long.boot.pp$Y_neg <- ifelse(dat.long.boot.pp$C==1, NA, dat.long.boot.pp$Y_neg)

# fit stage 1
prox_pooled_pp_s1 <- speedglm(Y_neg ~ ns(time_end, knots = c(10,20,30,40,50))*(treatment +
                                                                            # demographic
                                                                            sex_admin + age_years + bmi + race + charlson_cat_fac +
                                                                            # other
                                                                            ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks + 
                                                                            # NEC
                                                                            flu_vax),
                               data=dat.long.boot.pp,
                               family=binomial(),
                              sparse = FALSE)

# generate predictions from stage 1
dat.long.boot.pp$p_pp <- predict(prox_pooled_pp_s1, newdata = dat.long.boot.pp)

# fit stage 2
prox_pooled_pp_s2 <- speedglm(Y_pos ~ ns(time_end, knots = c(10,20,30,40,50))*treatment +
                                 # demographic
                                 sex_admin + age_years + bmi + race + charlson_cat_fac +
                                 # other
                                 ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks + 
                                 # predictions from stage 1
                                 p_pp,
                               # no NEC
                               data=dat.long.boot.pp,
                               family=binomial(),
                              sparse = FALSE)

# fit observed data model
prox_pooled_pp_obs <- speedglm(Y_pos ~ ns(time_end, knots = c(10,20,30,40,50))*treatment +
                                  # demographic
                                  sex_admin + age_years + bmi + race + charlson_cat_fac +
                                  # other
                                  ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks + 
                                  # NEC
                                  flu_vax,
                                data=dat.long.boot.pp,
                                family=binomial(),
                               sparse = FALSE)

# g formula setup

dat.boot$gmaxt <- 53

prox_pp_A0.long <- dat.boot[rep(1:nrow(dat.boot), dat.boot$gmaxt),]
prox_pp_A0.long$time_start <- ave(prox_pp_A0.long$bootid_mrn, prox_pp_A0.long$bootid_mrn, FUN=seq_along)
prox_pp_A0.long$time_start <- (prox_pp_A0.long$time_start-1)*time_unit
prox_pp_A0.long$time_end <- prox_pp_A0.long$time_start+time_unit
prox_pp_A0.long$treatment_obs <- prox_pp_A0.long$treatment
prox_pp_A0.long$treatment <- 0

prox_pp_A1.long <- dat.boot[rep(1:nrow(dat.boot), dat.boot$gmaxt),]
prox_pp_A1.long$time_start <- ave(prox_pp_A1.long$bootid_mrn, prox_pp_A1.long$bootid_mrn, FUN=seq_along)
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
prox_pp_A0.long$hazard_pos <- prox_pp_A0.long$hazard_pos_obs * exp(-prox_pp_A0.long$logHR * prox_pp_A0.long$treatment_obs)
### adding treated to untreated
prox_pp_A1.long$hazard_pos <- prox_pp_A1.long$hazard_pos_obs * exp(prox_pp_A1.long$logHR * (1-prox_pp_A1.long$treatment_obs))

### compute survival and cumulative incidence

### calculate (1 - hazard POSITIVE)
prox_pp_A0.long$pnoevent_pos <- 1 - prox_pp_A0.long$hazard_pos
prox_pp_A1.long$pnoevent_pos <- 1 - prox_pp_A1.long$hazard_pos

### calculate (1 - hazard NEGATIVE)
prox_pp_A0.long$pnoevent_neg <- 1 - prox_pp_A0.long$hazard_neg
prox_pp_A1.long$pnoevent_neg <- 1 - prox_pp_A1.long$hazard_neg

### sort the data by ID, time
prox_pp_A0.long <- prox_pp_A0.long[order(prox_pp_A0.long$bootid_mrn, prox_pp_A0.long$time_end),] 
prox_pp_A1.long <- prox_pp_A1.long[order(prox_pp_A1.long$bootid_mrn, prox_pp_A1.long$time_end),]

### lag (1 - hazard POSITIVE)
prox_pp_A0.long <- prox_pp_A0.long |> 
  arrange(bootid_mrn, time_end) |> 
  group_by(bootid_mrn) |> 
  mutate(pnoevent_pos_lag = lag(pnoevent_pos, n=1, default=1)) |> 
  ungroup()

prox_pp_A1.long <- prox_pp_A1.long |> 
  arrange(bootid_mrn, time_end) |> 
  group_by(bootid_mrn) |> 
  mutate(pnoevent_pos_lag = lag(pnoevent_pos, n=1, default=1)) |> 
  ungroup()

### lag (1 - hazard NEGATIVE)
prox_pp_A0.long <- prox_pp_A0.long |> 
  arrange(bootid_mrn, time_end) |> 
  group_by(bootid_mrn) |> 
  mutate(pnoevent_neg_lag = lag(pnoevent_neg, n=1, default=1)) |> 
  ungroup()

prox_pp_A1.long <- prox_pp_A1.long |> 
  arrange(bootid_mrn, time_end) |> 
  group_by(bootid_mrn) |> 
  mutate(pnoevent_neg_lag = lag(pnoevent_neg, n=1, default=1)) |> 
  ungroup()

### Aalen–Johansen estimator

prox_pp_A0.long$surv_prod <- prox_pp_A0.long$pnoevent_neg * prox_pp_A0.long$pnoevent_pos_lag 
prox_pp_A1.long$surv_prod <- prox_pp_A1.long$pnoevent_neg * prox_pp_A1.long$pnoevent_pos_lag

prox_pp_A0.long$survival <- ave(prox_pp_A0.long$surv_prod, prox_pp_A0.long$bootid_mrn, FUN=cumprod)
prox_pp_A1.long$survival <- ave(prox_pp_A1.long$surv_prod, prox_pp_A1.long$bootid_mrn, FUN=cumprod)

prox_pp_A0.long$risk_prod <- prox_pp_A0.long$hazard_pos * prox_pp_A0.long$survival
prox_pp_A1.long$risk_prod <- prox_pp_A1.long$hazard_pos * prox_pp_A1.long$survival

prox_pp_A0.long$risk_pos <- ave(prox_pp_A0.long$risk_prod, prox_pp_A0.long$bootid_mrn, FUN=cumsum)
prox_pp_A1.long$risk_pos <- ave(prox_pp_A1.long$risk_prod, prox_pp_A1.long$bootid_mrn, FUN=cumsum)

prox_pp_A0.long.res <- aggregate(risk_pos ~ time_end, data=prox_pp_A0.long, FUN=mean)
prox_pp_A1.long.res <- aggregate(risk_pos ~ time_end, data=prox_pp_A1.long, FUN=mean)

elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
message("Finished bootstrap ", i, " in ", round(elapsed, 2), " minutes")

saveRDS(
  cbind(sim = i,
        time_end = prox_pp_A0.long.res$time_end,
        risk0 = prox_pp_A0.long.res$risk_pos,
        risk1 = prox_pp_A1.long.res$risk_pos),
  file.path(res_path, sprintf("pci_pp_boot_rep_%03d.rds", i))
)

