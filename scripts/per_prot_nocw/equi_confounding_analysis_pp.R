##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Equi Confounding Analysis
##----- Per-protocol, no censoring weights

# Packages ----------------------------------------------------------------


library(tidyverse)
library(tidycmprsk)
library(survival)
library(ggsurvfit)
library(riskRegression)
library(geepack)
library(data.table)


# Data --------------------------------------------------------------------


data_Y3 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch/data_Y3_weekmatch.rds")

data_Y3 <- data_Y3 |> 
  mutate(Y3_pp_factor = case_when(
    Y3_pp_trunc==0 ~ "Censor",
    Y3_pp_trunc==1 ~ "Test Negative",
    Y3_pp_trunc==2 ~ "Test Positive"
  )) |> 
  mutate(Y3_pp_factor = factor(Y3_pp_factor, levels = c("Censor", "Test Negative", "Test Positive")),
         subclass=as.character(subclass))

res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results_weekmatch_pp/"


# Cox Model ---------------------------------------------------------------


# Per-protocol, no censoring weights

eqc_pp_fit1 <- coxph(
  Surv(Y3_pp_t_trunc, Y3_pp_factor == "Test Negative") ~ treatment +
    # demog
    sex_admin + age_years + bmi + race + charlson_cat_fac +
    # other
    ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
    # NEC
    flu_vax + 
    cluster(subclass),
  data = data_Y3
)


eqc_pp_fit2 <- coxph(
  Surv(Y3_pp_t_trunc, Y3_pp_factor == "Test Positive") ~ treatment +
    # demog
    sex_admin + age_years + bmi + race + charlson_cat_fac +
    # other
    ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
    # NEC
    flu_vax + 
    cluster(subclass),
  data = data_Y3
)

eqc.pp.cox.pointest <- c(
  treatHR= unname( exp(eqc_pp_fit2$coefficients["treatment"]) / exp(eqc_pp_fit1$coefficients["treatment"]) ),
  fluvaxHR=unname( exp(eqc_pp_fit2$coefficients["flu_vax"]) / exp(eqc_pp_fit1$coefficients["flu_vax"]) )
)

saveRDS(eqc.pp.cox.pointest, paste0(res_path,"eqc.pp.cox.pointest.rds")) # 2026-07-01


# Bootstrap CIs for PP ---------------------------------------------------


num.boot <- 100

set.seed(1155)
seed <- floor(runif(num.boot)*10^8)

setkey(data_Y3, subclass)
subclasses <- data_Y3[, unique(subclass)]
n_sub <- length(subclasses)


boot.results <- lapply(1:num.boot, function(i){
  
  t0 <- Sys.time()
  
  set.seed(seed[i])
  
  message("Starting EQC PP Cox bootstrap ", i, " (seed=", seed[i], ")")
  
  # select matched pairs
  samp_sub <- sample(subclasses, size = n_sub, replace = TRUE)
  
  # build boot dataset efficiently via one join:
  # map draw index j -> sampled subclass, then join to replicate all rows per subclass
  map <- data.table(j = seq_along(samp_sub), subclass = samp_sub)
  dat.boot <- data_Y3[map, on = "subclass", allow.cartesian = TRUE]
  # new cluster id per draw
  dat.boot[, subclass_boot := paste0(subclass, ".", j)] 
  
  # run test neg model
  eqc_pp_fit1 <- coxph(
    Surv(Y3_pp_t_trunc, Y3_pp_factor == "Test Negative") ~ treatment +
      # demog
      sex_admin + age_years + bmi + race + charlson_cat_fac +
      # other
      ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
      # NEC
      flu_vax + 
      cluster(subclass_boot),
    data = dat.boot
  )
  
  # run test pos model
  eqc_pp_fit2 <- coxph(
    Surv(Y3_pp_t_trunc, Y3_pp_factor == "Test Positive") ~ treatment +
      # demog
      sex_admin + age_years + bmi + race + charlson_cat_fac +
      # other
      ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
      # NEC
      flu_vax + 
      cluster(subclass_boot),
    data = dat.boot
  )
  
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  message("Finished bootstrap ", i, " in ", round(elapsed, 2), " minutes")
  
  return(c(
    sim=i,
    treatHR=unname(exp(eqc_pp_fit2$coefficients["treatment"])) / unname(exp(eqc_pp_fit1$coefficients["treatment"])),
    fluvaxHR=unname(exp(eqc_pp_fit2$coefficients["flu_vax"])) / unname(exp(eqc_pp_fit1$coefficients["flu_vax"]))
  ))
  
})

boot.long <- bind_rows(lapply(boot.results, tibble::as_tibble_row))
saveRDS(boot.long, paste0(res_path, "eqc.pp.cox.boot.long.rds"))




