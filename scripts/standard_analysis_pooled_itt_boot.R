##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Standard Analysis Pooled ITT Bootstrap

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(data.table)
library(speedglm)
library(splines)
library(geepack)

# Data --------------------------------------------------------------------

data_Y2 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_Y2.rds")
dat <- data_Y2

res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results/"

# Boot --------------------------------------------------------------------

num.boot <- 100

set.seed(1155)
seed <- floor(runif(num.boot)*10^8)

boot.results <- lapply(1:num.boot, function(i){
  
  t0 <- Sys.time()
  
  set.seed(seed[i])
  
  message("Starting bootstrap ", i, " (seed=", seed[i], ")")
  
  # select row indices
  sample.row <- sample(nrow(dat), size=nrow(dat), replace=T)
  
  # pull selected rows
  dat.boot <- dat[sample.row,]  
  dat.boot$bootid <- seq(1:nrow(dat.boot))
  
  # long format data
  time_unit <- 1
  
  dat.boot$max_units <- ceiling(dat.boot$Y2_itt_t_trunc/time_unit)+1
  dat.long.boot.itt <- dat.boot[rep(1:nrow(dat.boot), dat.boot$max_units),]
  
  dat.long.boot.itt$time_start <- ave(dat.long.boot.itt$bootid, dat.long.boot.itt$bootid, FUN=seq_along)
  dat.long.boot.itt$time_start <- (dat.long.boot.itt$time_start-1)*time_unit
  dat.long.boot.itt$time_end <- dat.long.boot.itt$time_start+time_unit
  
  dat.long.boot.itt$Y <- ifelse(
    dat.long.boot.itt$Y2_itt_trunc == 1 &
      dat.long.boot.itt$Y2_itt_t_trunc == dat.long.boot.itt$time_start,
    1, 0
  )
  
  dat.long.boot.itt$C <- ifelse(
    dat.long.boot.itt$Y2_itt_trunc == 0 &
      dat.long.boot.itt$Y2_itt_t_trunc == dat.long.boot.itt$time_start,
    1, 0
  )
  
  dat.long.boot.itt$Y <- ifelse(dat.long.boot.itt$C==1, NA, dat.long.boot.itt$Y)
  
  # run regression
  std_pooled_itt <- speedglm(Y ~ ns(time_end, knots = c(10,20,30))*treatment +
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
  
  std_itt_A0.long <- dat.boot[rep(1:nrow(dat.boot), dat.boot$gmaxt),]
  std_itt_A0.long$time_start <- ave(std_itt_A0.long$bootid, std_itt_A0.long$bootid, FUN=seq_along)
  std_itt_A0.long$time_start <- (std_itt_A0.long$time_start-1)*time_unit
  std_itt_A0.long$time_end <- std_itt_A0.long$time_start+time_unit
  std_itt_A0.long$treatment <- 0
  
  std_itt_A1.long <- dat.boot[rep(1:nrow(dat.boot), dat.boot$gmaxt),]
  std_itt_A1.long$time_start <- ave(std_itt_A1.long$bootid, std_itt_A1.long$bootid, FUN=seq_along)
  std_itt_A1.long$time_start <- (std_itt_A1.long$time_start-1)*time_unit
  std_itt_A1.long$time_end <- std_itt_A1.long$time_start+time_unit
  std_itt_A1.long$treatment <- 1
  
  std_itt_A0.long$hazard <- predict(std_pooled_itt, newdata=std_itt_A0.long, type="response")
  std_itt_A1.long$hazard <- predict(std_pooled_itt, newdata=std_itt_A1.long, type="response")
  
  std_itt_A0.long$pnoevent <- 1 - std_itt_A0.long$hazard
  std_itt_A1.long$pnoevent <- 1 - std_itt_A1.long$hazard
  
  std_itt_A0.long <- std_itt_A0.long[order(std_itt_A0.long$bootid, std_itt_A0.long$time_end),]
  std_itt_A1.long <- std_itt_A1.long[order(std_itt_A1.long$bootid, std_itt_A1.long$time_end),]
  
  std_itt_A0.long$survival <- ave(std_itt_A0.long$pnoevent, std_itt_A0.long$bootid, FUN=cumprod)
  std_itt_A1.long$survival <- ave(std_itt_A1.long$pnoevent, std_itt_A1.long$bootid, FUN=cumprod)
  
  std_itt_A0.long$risk <- 1 - std_itt_A0.long$survival
  std_itt_A1.long$risk <- 1 - std_itt_A1.long$survival
  
  std_itt_A0.long <- aggregate(risk ~ time_end, data=std_itt_A0.long, FUN=mean)
  std_itt_A1.long <- aggregate(risk ~ time_end, data=std_itt_A1.long, FUN=mean)
  
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  message("Finished bootstrap ", i, " in ", round(elapsed, 2), " minutes")
  
  return(cbind(sim=i, 
               time_end=std_itt_A0.long$time_end,
               risk0=std_itt_A0.long$risk,
               risk1=std_itt_A1.long$risk))
  
})

boot.long <- bind_rows(lapply(boot.results, as.data.frame))
saveRDS(boot.long, paste0(res_path, "std.itt.boot.long.rds"))

