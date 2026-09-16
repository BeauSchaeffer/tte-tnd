##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Proximal Inference Analysis Pooled PP Bootstrap -- NET RISK, Option B
##----- Per-protocol, no censoring weights
##----- ** SINGLE BOOTSTRAP REPLICATE FOR USE WITH ARRAY **
##----- last updated 2026-09-16


# Packages ----------------------------------------------------------------


library(tidyverse)
library(data.table)
library(speedglm)
library(splines)


# Data --------------------------------------------------------------------


data_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch.4/"

data_Y2  <- read_rds(paste0(data_path, "data_Y2_weekmatch.rds"))
neg_hist <- read_rds(paste0(data_path, "neg_hist_weekmatch.rds"))

dat <- data_Y2
setDT(dat)

res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results_pp.5/netB_pci_boot_reps"
dir.create(res_path, showWarnings = FALSE, recursive = TRUE)

### Negative-test weeks inside each person's PP risk set. Does not depend on
### the resample; joined on (fake_mrn, time_start) below, so every drawn copy
### of a person receives the same negative weeks.

neg_wk <- neg_hist |>
  inner_join(dat |> dplyr::select(fake_mrn, cap = Y2_pp_t_trunc), by = "fake_mrn") |>
  filter(neg_t <= cap) |>
  distinct(fake_mrn, time_start = neg_t)

setDT(neg_wk)


# Boot --------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])

num.boot <- 200

stopifnot(!is.na(i), i >= 1, i <= num.boot)

set.seed(1155)
seed <- floor(runif(num.boot)*10^8)
set.seed(seed[i])

setDT(dat)
setkey(dat, subclass)
subclasses <- dat[, unique(subclass)]
n_sub <- length(subclasses)

t0 <- Sys.time()

message(
  "Starting PCI net-B PP pooled bootstrap ", i,
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


# PP long format expansion -----------------------------------------------


### calc number of rows needed for each individual
time_unit <- 1

### ensure at least 1 row for each individual
dat.boot$max_units <- ceiling(dat.boot$Y2_pp_t_trunc/time_unit)+1
dat.long.boot.pp <- dat.boot[rep(1:nrow(dat.boot), dat.boot$max_units),]

### variable that represents the start and end time corresponding to each row of observation
dat.long.boot.pp$time_start <- ave(dat.long.boot.pp$bootid_mrn, dat.long.boot.pp$bootid_mrn, FUN=seq_along)
dat.long.boot.pp$time_start <- (dat.long.boot.pp$time_start-1)*time_unit
dat.long.boot.pp$time_end <- dat.long.boot.pp$time_start+time_unit

# recommended add
dat.long.boot.pp <- dat.long.boot.pp[order(dat.long.boot.pp$bootid_mrn, dat.long.boot.pp$time_end),]

### modify the Y and C variables so that they are only equal to 1 if the 
### event/censoring happened in that time interval
dat.long.boot.pp$Y_pos <- ifelse(
  dat.long.boot.pp$Y2_pp_trunc == 1 &
    dat.long.boot.pp$Y2_pp_t_trunc == dat.long.boot.pp$time_start,
  1, 0
)

dat.long.boot.pp$C <- ifelse(
  dat.long.boot.pp$Y2_pp_trunc == 0 &
    dat.long.boot.pp$Y2_pp_t_trunc == dat.long.boot.pp$time_start,
  1, 0
)

### Y_neg is recurrent: 1 in every week a negative test occurred, marked on
### the same rows as Y_pos. Negatives at the cap are kept (<=) since time is
### in whole weeks; same-week repeats collapse to a single indicator.

setDT(dat.long.boot.pp)

dat.long.boot.pp[, Y_neg := 0]
dat.long.boot.pp[neg_wk, on = c("fake_mrn", "time_start"), Y_neg := 1]

dat.long.boot.pp$Y_pos <- ifelse(dat.long.boot.pp$C==1, NA, dat.long.boot.pp$Y_pos)
dat.long.boot.pp$Y_neg <- ifelse(dat.long.boot.pp$C==1, NA, dat.long.boot.pp$Y_neg)



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
                               data=dat.long.boot.pp,
                               family=binomial(),
                              sparse = FALSE)

dat.long.boot.pp$p_pp <- predict(prox_pooled_pp_s1, newdata = dat.long.boot.pp)

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


# PP Survival and Risk ---------------------------------------------------

### Net risk is a SINGLE decrement: only the positive hazard removes anyone
### from the risk set, so survival is prod(1 - hazard_pos) alone. The negative
### hazard does not enter the survival product -- under PCI it enters only
### through the stage 1 bridge p_pp.

dat.boot$gmaxt <- 53

### G formula data setup A=0
prox_pp_A0.long <- dat.boot[rep(1:nrow(dat.boot), dat.boot$gmaxt),]
prox_pp_A0.long$time_start <- ave(prox_pp_A0.long$bootid_mrn, prox_pp_A0.long$bootid_mrn, FUN=seq_along)
prox_pp_A0.long$time_start <- (prox_pp_A0.long$time_start-1)*time_unit
prox_pp_A0.long$time_end <- prox_pp_A0.long$time_start+time_unit
prox_pp_A0.long$treatment_obs <- prox_pp_A0.long$treatment
prox_pp_A0.long$treatment <- 0

### G formula data setup A=1
prox_pp_A1.long <- dat.boot[rep(1:nrow(dat.boot), dat.boot$gmaxt),]
prox_pp_A1.long$time_start <- ave(prox_pp_A1.long$bootid_mrn, prox_pp_A1.long$bootid_mrn, FUN=seq_along)
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


### Kaplan-Meier type estimator (single decrement -- lagged positive hazard
### only, no competing test-negative term)

prox_pp_A0.long$survival <- ave(prox_pp_A0.long$pnoevent_pos_lag, prox_pp_A0.long$bootid_mrn, FUN=cumprod)
prox_pp_A1.long$survival <- ave(prox_pp_A1.long$pnoevent_pos_lag, prox_pp_A1.long$bootid_mrn, FUN=cumprod)

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
  file.path(res_path, sprintf("netB_pci_pp_boot_rep_%03d.rds", i))
)
