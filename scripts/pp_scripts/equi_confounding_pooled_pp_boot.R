##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Equi Confounding Analysis Pooled PP Bootstrap
##----- Per-protocol, no censoring weights

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(data.table)
library(speedglm)
library(splines)


# Data --------------------------------------------------------------------


data_Y3 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch/data_Y3_weekmatch.rds")
dat <- data_Y3
setDT(dat)

res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results_weekmatch_pp/"


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
    "Starting EQC PP pooled bootstrap ", i,
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
  eqc_pooled_pp_fit1 <- speedglm(Y_neg ~ ns(time_end, knots = c(10,20,30))*treatment +
                                    # demographic
                                    sex_admin + age_years + bmi + race + charlson_cat_fac +
                                    # other
                                    ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks + 
                                    # NEC
                                    flu_vax,
                                  data=dat.long.boot.pp,
                                  family=binomial())
  
  # fit stage 2
  eqc_pooled_pp_fit2 <- speedglm(Y_pos ~ ns(time_end, knots = c(10,20,30))*treatment +
                                    # demographic
                                    sex_admin + age_years + bmi + race + charlson_cat_fac +
                                    # other
                                    ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks + 
                                    # NEC
                                    flu_vax,
                                  data=dat.long.boot.pp,
                                  family=binomial())
  
  # g formula setup
  
  dat.boot$gmaxt <- 53
  
  eqc_pp_A0.long <- dat.boot[rep(1:nrow(dat.boot), dat.boot$gmaxt),]
  eqc_pp_A0.long$time_start <- ave(eqc_pp_A0.long$bootid_mrn, eqc_pp_A0.long$bootid_mrn, FUN=seq_along)
  eqc_pp_A0.long$time_start <- (eqc_pp_A0.long$time_start-1)*time_unit
  eqc_pp_A0.long$time_end <- eqc_pp_A0.long$time_start+time_unit
  eqc_pp_A0.long$treatment <- 0
  
  eqc_pp_A1.long <- dat.boot[rep(1:nrow(dat.boot), dat.boot$gmaxt),]
  eqc_pp_A1.long$time_start <- ave(eqc_pp_A1.long$bootid_mrn, eqc_pp_A1.long$bootid_mrn, FUN=seq_along)
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
  eqc_pp_A0.long <- eqc_pp_A0.long[order(eqc_pp_A0.long$bootid_mrn, eqc_pp_A0.long$time_end),] 
  eqc_pp_A1.long <- eqc_pp_A1.long[order(eqc_pp_A1.long$bootid_mrn, eqc_pp_A1.long$time_end),] 
  
  ### Calculate the cumulative survival 
  
  # lag P(no event pos)
  
  eqc_pp_A0.long <- eqc_pp_A0.long |> 
    arrange(bootid_mrn, time_end) |> 
    group_by(bootid_mrn) |> 
    mutate(pnoevent_pos_lag = lag(pnoevent_pos, n=1, default=1),
           pnoevent_pos_c_lag = lag(pnoevent_pos_c, n=1, default=1)) |> 
    ungroup()
  
  eqc_pp_A1.long <- eqc_pp_A1.long |> 
    arrange(bootid_mrn, time_end) |> 
    group_by(bootid_mrn) |> 
    mutate(pnoevent_pos_lag = lag(pnoevent_pos, n=1, default=1)) |> 
    ungroup()
  
  # product at each time (P(no event neg) * lag P(no event pos))
  
  eqc_pp_A0.long$surv_prod_lag <- eqc_pp_A0.long$pnoevent_neg * eqc_pp_A0.long$pnoevent_pos_lag
  eqc_pp_A1.long$surv_prod_lag <- eqc_pp_A1.long$pnoevent_neg * eqc_pp_A1.long$pnoevent_pos_lag
  eqc_pp_A0.long$surv_prod_c_lag <- eqc_pp_A0.long$pnoevent_neg * eqc_pp_A0.long$pnoevent_pos_c_lag
  
  # cumulative product within individual
  
  eqc_pp_A0.long$survival_pos <- ave(eqc_pp_A0.long$surv_prod_lag, eqc_pp_A0.long$bootid_mrn, FUN=cumprod)
  eqc_pp_A1.long$survival_pos <- ave(eqc_pp_A1.long$surv_prod_lag, eqc_pp_A1.long$bootid_mrn, FUN=cumprod)
  eqc_pp_A0.long$survival_pos_c <- ave(eqc_pp_A0.long$surv_prod_c_lag, eqc_pp_A0.long$bootid_mrn, FUN=cumprod)
  
  ### Calculate risk using CIF estimator
  
  # product at each time (haz pos * surv pos)
  
  eqc_pp_A0.long$risk_prod_pos <- eqc_pp_A0.long$hazard_pos * eqc_pp_A0.long$survival_pos
  eqc_pp_A1.long$risk_prod_pos <- eqc_pp_A1.long$hazard_pos * eqc_pp_A1.long$survival_pos
  eqc_pp_A0.long$risk_prod_pos_c <- eqc_pp_A0.long$hazard_pos_c * eqc_pp_A0.long$survival_pos_c
  
  # cumulative sum within individual
  
  eqc_pp_A0.long$risk_pos <- ave(eqc_pp_A0.long$risk_prod_pos, eqc_pp_A0.long$bootid_mrn, FUN=cumsum)
  eqc_pp_A1.long$risk_pos <- ave(eqc_pp_A1.long$risk_prod_pos, eqc_pp_A1.long$bootid_mrn, FUN=cumsum)
  eqc_pp_A0.long$risk_pos_c <- ave(eqc_pp_A0.long$risk_prod_pos_c, eqc_pp_A0.long$bootid_mrn, FUN=cumsum)
  
  # Calculate the average risk at each time point
  
  eqc_pp_A0.long.res <- aggregate(risk_pos ~ time_end, data=eqc_pp_A0.long, FUN=mean)
  eqc_pp_A1.long.res <- aggregate(risk_pos ~ time_end, data=eqc_pp_A1.long, FUN=mean)
  eqc_pp_A0.long.res.c <- aggregate(risk_pos_c ~ time_end, data=eqc_pp_A0.long, FUN=mean)
  
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  message("Finished bootstrap ", i, " in ", round(elapsed, 2), " minutes")
  
  return(cbind(sim = i, 
               time_end = eqc_pp_A0.long.res$time_end,
               risk0 = eqc_pp_A0.long.res$risk_pos,
               risk0corr = eqc_pp_A0.long.res.c$risk_pos_c,
               risk1 = eqc_pp_A1.long.res$risk_pos))
})

boot.long <- bind_rows(lapply(boot.results, as.data.frame))
saveRDS(boot.long, paste0(res_path, "eqc.pp.boot.long.rds"))














