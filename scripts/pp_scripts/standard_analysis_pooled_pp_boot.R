##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Standard Analysis Pooled PP Bootstrap
##----- Per-protocol, no censoring weights
##----- last updated 2026-07-16

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(data.table)
library(speedglm)
library(splines)

# Data --------------------------------------------------------------------

data_Y2 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch/data_Y2_weekmatch.rds")
dat <- data_Y2

res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results_pp.2/"

# Boot --------------------------------------------------------------------

num.boot <- 50

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
    "Starting STD PP pooled bootstrap ", i,
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
  
  dat.boot$max_units <- ceiling(dat.boot$Y2_pp_t_trunc/time_unit)+1
  dat.long.boot.pp <- dat.boot[rep(1:nrow(dat.boot), dat.boot$max_units),]
  
  dat.long.boot.pp$time_start <- ave(dat.long.boot.pp$bootid_mrn, dat.long.boot.pp$bootid_mrn, FUN=seq_along)
  dat.long.boot.pp$time_start <- (dat.long.boot.pp$time_start-1)*time_unit
  dat.long.boot.pp$time_end <- dat.long.boot.pp$time_start+time_unit
  
  # recommended add
  # dat.long.boot.pp <- dat.long.boot.pp[order(dat.long.boot.pp$bootid_mrn, dat.long.boot.pp$time_end),]
  
  dat.long.boot.pp$Y <- ifelse(
    dat.long.boot.pp$Y2_pp_trunc == 1 &
      dat.long.boot.pp$Y2_pp_t_trunc == dat.long.boot.pp$time_start,
    1, 0
  )
  
  dat.long.boot.pp$C <- ifelse(
    dat.long.boot.pp$Y2_pp_trunc == 0 &
      dat.long.boot.pp$Y2_pp_t_trunc == dat.long.boot.pp$time_start,
    1, 0
  )
  
  dat.long.boot.pp$Y <- ifelse(dat.long.boot.pp$C==1, NA, dat.long.boot.pp$Y)
  
  # run regression
  std_pooled_pp <- speedglm(Y ~ ns(time_end, knots = c(10,20,30,40,50))*treatment +
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
  
  std_pp_A0.long <- dat.boot[rep(1:nrow(dat.boot), dat.boot$gmaxt),]
  std_pp_A0.long$time_start <- ave(std_pp_A0.long$bootid_mrn, std_pp_A0.long$bootid_mrn, FUN=seq_along)
  std_pp_A0.long$time_start <- (std_pp_A0.long$time_start-1)*time_unit
  std_pp_A0.long$time_end <- std_pp_A0.long$time_start+time_unit
  std_pp_A0.long$treatment <- 0
  
  std_pp_A1.long <- dat.boot[rep(1:nrow(dat.boot), dat.boot$gmaxt),]
  std_pp_A1.long$time_start <- ave(std_pp_A1.long$bootid_mrn, std_pp_A1.long$bootid_mrn, FUN=seq_along)
  std_pp_A1.long$time_start <- (std_pp_A1.long$time_start-1)*time_unit
  std_pp_A1.long$time_end <- std_pp_A1.long$time_start+time_unit
  std_pp_A1.long$treatment <- 1
  
  std_pp_A0.long$hazard <- predict(std_pooled_pp, newdata=std_pp_A0.long, type="response")
  std_pp_A1.long$hazard <- predict(std_pooled_pp, newdata=std_pp_A1.long, type="response")
  
  std_pp_A0.long$pnoevent <- 1 - std_pp_A0.long$hazard
  std_pp_A1.long$pnoevent <- 1 - std_pp_A1.long$hazard
  
  std_pp_A0.long <- std_pp_A0.long[order(std_pp_A0.long$bootid_mrn, std_pp_A0.long$time_end),]
  std_pp_A1.long <- std_pp_A1.long[order(std_pp_A1.long$bootid_mrn, std_pp_A1.long$time_end),]
  
  std_pp_A0.long$survival <- ave(std_pp_A0.long$pnoevent, std_pp_A0.long$bootid_mrn, FUN=cumprod)
  std_pp_A1.long$survival <- ave(std_pp_A1.long$pnoevent, std_pp_A1.long$bootid_mrn, FUN=cumprod)
  
  std_pp_A0.long$risk <- 1 - std_pp_A0.long$survival
  std_pp_A1.long$risk <- 1 - std_pp_A1.long$survival
  
  std_pp_A0.long <- aggregate(risk ~ time_end, data=std_pp_A0.long, FUN=mean)
  std_pp_A1.long <- aggregate(risk ~ time_end, data=std_pp_A1.long, FUN=mean)
  
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  message("Finished bootstrap ", i, " in ", round(elapsed, 2), " minutes")
  
  return(cbind(sim=i, 
               time_end=std_pp_A0.long$time_end,
               risk0=std_pp_A0.long$risk,
               risk1=std_pp_A1.long$risk))
  
})

boot.long <- bind_rows(lapply(boot.results, as.data.frame))
saveRDS(boot.long, paste0(res_path, "std.pp.boot.long.rds"))

