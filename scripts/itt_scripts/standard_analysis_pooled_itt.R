##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Standard Analysis Pooled
##----- last updated 2026-07-09

# Packages ----------------------------------------------------------------


library(tidyverse)
library(data.table)
library(speedglm)
library(splines)
library(geepack)


# Data --------------------------------------------------------------------


# data_Y2 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_Y2.rds")
data_Y2 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch/data_Y2_weekmatch.rds")
dat <- data_Y2

res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results_itt.2/"

# Downsample --------------------------------------------------------------


# subclass_ids <- data_Y2 |> dplyr::select(subclass) |> unique()
# set.seed(345)
# subclass_ids_subset <- slice_sample(subclass_ids, n=150000) # works with 150000
# dat_downsamp <- data_Y2 |> filter(subclass %in% subclass_ids_subset$subclass) |> droplevels()
# rm(subclass_ids, subclass_ids_subset)


# ITT long format expansion -----------------------------------------------


  ### calc number of rows needed for each individual
time_unit <- 1
  
### ensure at least 1 row for each individual
dat$max_units <- ceiling(dat$Y2_itt_t_trunc/time_unit)+1
dat.long.itt <- dat[rep(1:nrow(dat), dat$max_units),]
  
### variable that represents the start and end time corresponding to each row of observation
dat.long.itt$time_start <- ave(dat.long.itt$fake_mrn, dat.long.itt$fake_mrn, FUN=seq_along)
dat.long.itt$time_start <- (dat.long.itt$time_start-1)*time_unit
dat.long.itt$time_end <- dat.long.itt$time_start+time_unit
  
### modify the Y and C variables so that they are only equal to 1 if the 
  ### event/censoring happened in that time interval
dat.long.itt$Y <- ifelse(
  dat.long.itt$Y2_itt_trunc == 1 &
    dat.long.itt$Y2_itt_t_trunc == dat.long.itt$time_start,
  1, 0
)

dat.long.itt$C <- ifelse(
  dat.long.itt$Y2_itt_trunc == 0 &
    dat.long.itt$Y2_itt_t_trunc == dat.long.itt$time_start,
  1, 0
)

dat.long.itt$Y <- ifelse(dat.long.itt$C==1, NA, dat.long.itt$Y)


# ITT Pooled Logistic -----------------------------------------------------


std_pooled_itt <- speedglm(Y ~ ns(time_end, knots = c(10,20,30,40,50))*treatment +
                             # demographic
                             sex_admin + age_years + bmi + race + charlson_cat_fac +
                             # other
                             ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks + 
                             # NEC
                             flu_vax,
                         data=dat.long.itt,
                         family=binomial())

summary(std_pooled_itt)
saveRDS(std_pooled_itt, paste0(res_path,"std_pooled_itt_model.rds"))

  ### sanity check
  ### fit same model without interaction terms
  ### coefficients should be equivalent/similar to Cox

# std_pooled_itt_noint <- speedglm(Y ~ ns(time_end, knots = c(10,20,30,40,50)) + treatment +
#                              # demographic
#                              sex_admin + age_years + bmi + race + charlson_cat_fac +
#                              # other
#                              ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks + 
#                              # NEC
#                              flu_vax,
#                            data=dat.long.itt,
#                            family=binomial())
# 
# std_pooled_itt_noint_tidy <- geepack::tidy(std_pooled_itt_noint, conf.int = TRUE, exponentiate = TRUE)
# std_pooled_itt_noint_tidy |> print(n=100)


# ITT Survival and Risk ---------------------------------------------------


dat$gmaxt <- 53

  ### G formula data setup A=0
std_itt_A0.long <- dat[rep(1:nrow(dat), dat$gmaxt),]
std_itt_A0.long$time_start <- ave(std_itt_A0.long$fake_mrn, std_itt_A0.long$fake_mrn, FUN=seq_along)
std_itt_A0.long$time_start <- (std_itt_A0.long$time_start-1)*time_unit
std_itt_A0.long$time_end <- std_itt_A0.long$time_start+time_unit
std_itt_A0.long$treatment <- 0

  ### G formula data setup A=1
std_itt_A1.long <- dat[rep(1:nrow(dat), dat$gmaxt),]
std_itt_A1.long$time_start <- ave(std_itt_A1.long$fake_mrn, std_itt_A1.long$fake_mrn, FUN=seq_along)
std_itt_A1.long$time_start <- (std_itt_A1.long$time_start-1)*time_unit
std_itt_A1.long$time_end <- std_itt_A1.long$time_start+time_unit
std_itt_A1.long$treatment <- 1

### Calculate predicted hazards:
std_itt_A0.long$hazard <- predict(std_pooled_itt, newdata=std_itt_A0.long, type="response")
std_itt_A1.long$hazard <- predict(std_pooled_itt, newdata=std_itt_A1.long, type="response")

### Calculate (1 - hazard)
std_itt_A0.long$pnoevent <- 1 - std_itt_A0.long$hazard
std_itt_A1.long$pnoevent <- 1 - std_itt_A1.long$hazard

### Sort the data by time

  #* sort by ID, time

std_itt_A0.long <- std_itt_A0.long[order(std_itt_A0.long$time_end),]
std_itt_A1.long <- std_itt_A1.long[order(std_itt_A1.long$time_end),]

  # ### Calculate the cumulative survival by taking the cumulative product
  # ### of (1 - hazard)
  # std_itt_A0$survival <- cumprod(std_itt_A0$pnoevent)
  # std_itt_A1$survival <- cumprod(std_itt_A1$pnoevent)

std_itt_A0.long$survival <- ave(std_itt_A0.long$pnoevent, std_itt_A0.long$fake_mrn, FUN=cumprod)
std_itt_A1.long$survival <- ave(std_itt_A1.long$pnoevent, std_itt_A1.long$fake_mrn, FUN=cumprod)

### Calculate risk = 1 - survival
std_itt_A0.long$risk <- 1 - std_itt_A0.long$survival
std_itt_A1.long$risk <- 1 - std_itt_A1.long$survival

# Calculate the average risk at each time point. Here, we're going to use
# the aggregate function to do so:
std_itt_A0.long <- aggregate(risk ~ time_end, data=std_itt_A0.long, FUN=mean)
std_itt_A1.long <- aggregate(risk ~ time_end, data=std_itt_A1.long, FUN=mean)

# save point estimate risk curves in bootstrap-compatible format
std.itt.risk.pointest <- tibble(
  sim = 0L,  # 0 = main analysis (bootstraps are 1..B)
  time_end = std_itt_A0.long$time_end,
  risk0 = std_itt_A0.long$risk,
  risk1 = std_itt_A1.long$risk
)

saveRDS(std.itt.risk.pointest, paste0(res_path, "std.itt.risk.pointest.rds"))
