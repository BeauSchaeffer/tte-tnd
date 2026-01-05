##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Equi Confounding Analysis Pooled ITT Bootstrap

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(data.table)
library(speedglm)
library(splines)


# Data --------------------------------------------------------------------


data_Y3 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_Y3.rds")
dat <- data_Y3
setDT(dat)

res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results/"


# Boot --------------------------------------------------------------------

num.boot <- 100

set.seed(1155)
seed <- floor(runif(num.boot)*10^8)

setDT(dat)
setkey(dat, subclass)
subclasses <- dat[, unique(subclass)]
n_sub <- length(subclasses)

boot.results <- lapply(1:num.boot, function(i){
  
  t0 <- Sys.time()
  
  set.seed(seed[i])
  
  message(
    "Starting EQC ITT pooled bootstrap ", i,
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
  
  dat.boot$max_units <- ceiling(dat.boot$Y3_itt_t_trunc/time_unit)+1
  dat.long.boot.itt <- dat.boot[rep(1:nrow(dat.boot), dat.boot$max_units),]
  
  dat.long.boot.itt$time_start <- ave(dat.long.boot.itt$fake_mrn, dat.long.boot.itt$fake_mrn, FUN=seq_along)
  dat.long.boot.itt$time_start <- (dat.long.boot.itt$time_start-1)*time_unit
  dat.long.boot.itt$time_end <- dat.long.boot.itt$time_start+time_unit
  
  # recommended add
  dat.long.boot.itt <- dat.long.boot.itt[order(dat.long.boot.itt$bootid_mrn, dat.long.boot.itt$time_end),]
  
  dat.long.boot.itt$Y_pos <- ifelse(
    dat.long.boot.itt$Y3_itt_trunc == 2 &
      dat.long.boot.itt$Y3_itt_t_trunc == dat.long.boot.itt$time_start,
    1, 0
  )
  
  dat.long.boot.itt$Y_neg <- ifelse(
    dat.long.boot.itt$Y3_itt_trunc == 1 &
      dat.long.boot.itt$Y3_itt_t_trunc == dat.long.boot.itt$time_start,
    1, 0
  )
  
  dat.long.boot.itt$C <- ifelse(
    dat.long.boot.itt$Y3_itt_trunc == 0 &
      dat.long.boot.itt$Y3_itt_t_trunc == dat.long.boot.itt$time_start,
    1, 0
  )
  
  dat.long.boot.itt$Y_pos <- ifelse(dat.long.boot.itt$C==1, NA, dat.long.boot.itt$Y_pos)
  dat.long.boot.itt$Y_neg <- ifelse(dat.long.boot.itt$C==1, NA, dat.long.boot.itt$Y_neg)
  
  # fit stage 1
  eqc_pooled_itt_fit1 <- speedglm(Y_neg ~ ns(time_end, knots = c(10,20,30))*treatment +
                                    # demographic
                                    sex_admin + age_years + bmi + race + charlson_cat_fac +
                                    # other
                                    ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks + 
                                    # NEC
                                    flu_vax,
                                  data=dat.long.boot.itt,
                                  family=binomial())
  
  # fit stage 2
  eqc_pooled_itt_fit2 <- speedglm(Y_pos ~ ns(time_end, knots = c(10,20,30))*treatment +
                                    # demographic
                                    sex_admin + age_years + bmi + race + charlson_cat_fac +
                                    # other
                                    ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks + 
                                    # NEC
                                    flu_vax,
                                  data=dat.long.boot.itt,
                                  family=binomial())
  
  # g formula setup
  
  dat.boot$gmaxt <- 53
  
  ### G formula data setup A=0
  eqc_itt_A0.long <- dat.boot[rep(1:nrow(dat.boot), dat.boot$gmaxt),]
  eqc_itt_A0.long$time_start <- ave(eqc_itt_A0.long$fake_mrn, eqc_itt_A0.long$fake_mrn, FUN=seq_along)
  eqc_itt_A0.long$time_start <- (eqc_itt_A0.long$time_start-1)*time_unit
  eqc_itt_A0.long$time_end <- eqc_itt_A0.long$time_start+time_unit
  eqc_itt_A0.long$treatment <- 0
  
  ### G formula data setup A=1
  eqc_itt_A1.long <- dat.boot[rep(1:nrow(dat.boot), dat.boot$gmaxt),]
  eqc_itt_A1.long$time_start <- ave(eqc_itt_A1.long$fake_mrn, eqc_itt_A1.long$fake_mrn, FUN=seq_along)
  eqc_itt_A1.long$time_start <- (eqc_itt_A1.long$time_start-1)*time_unit
  eqc_itt_A1.long$time_end <- eqc_itt_A1.long$time_start+time_unit
  eqc_itt_A1.long$treatment <- 1
  
  # in progress
  
})















