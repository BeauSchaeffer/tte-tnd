##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Equi Confounding Analysis

# Packages ----------------------------------------------------------------


library(tidyverse)
library(tidycmprsk)
library(survival)
library(ggsurvfit)
library(riskRegression)
library(geepack)


# Data --------------------------------------------------------------------


data_Y3 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_Y3.rds")

data_Y3 <- data_Y3 |> 
  mutate(Y3_itt_factor = case_when(
    Y3_itt_trunc==0 ~ "Censor",
    Y3_itt_trunc==1 ~ "Test Negative",
    Y3_itt_trunc==2 ~ "Test Positive"
  )) |> 
  mutate(Y3_itt_factor = factor(Y3_itt_factor, levels = c("Censor", "Test Negative", "Test Positive")),
         subclass=as.character(subclass))

res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results/"


# Cox Model ---------------------------------------------------------------


# Intention to Treat

eqc_itt_fit1 <- coxph(
  Surv(Y3_itt_t_trunc, Y3_itt_factor == "Test Negative") ~ treatment +
    # demog
    sex_admin + age_years + bmi + race + charlson_cat_fac +
    # other
    ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
    # NEC
    flu_vax + 
    cluster(subclass),
  data = data_Y3
)


eqc_itt_fit2 <- coxph(
  Surv(Y3_itt_t_trunc, Y3_itt_factor == "Test Positive") ~ treatment +
    # demog
    sex_admin + age_years + bmi + race + charlson_cat_fac +
    # other
    ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
    # NEC
    flu_vax + 
    cluster(subclass),
  data = data_Y3
)

eqc.itt.cox.pointest <- c(
  exp(eqc_itt_fit2$coefficients[1]) / exp(eqc_itt_fit1$coefficients[1]),
  exp(eqc_itt_fit2$coefficients[21]) / exp(eqc_itt_fit1$coefficients[21])
  )
# saveRDS(eqc.itt.cox.pointest, paste0(res_path,"eqc.itt.cox.pointest.rds")) # 2025-12-23


# Bootstrap CIs for ITT ---------------------------------------------------


num.boot <- 2

set.seed(1155)
seed <- floor(runif(num.boot)*10^8)

by_sub <- split(data_Y3, data_Y3$subclass)
subclasses <- names(by_sub)
n_sub <- length(subclasses)


boot.results <- lapply(1:num.boot, function(i){
  
  t0 <- Sys.time()
  
  set.seed(seed[i])
  
  message("Starting EQC ITT Cox bootstrap ", i, " (seed=", seed[i], ")")
  
  # select matched pairs
  samp_sub <- sample(subclasses, size = n_sub, replace = TRUE)
  
  # build boot dataset: include all rows from each sampled subclass
  # and assign a fresh subclass ID so duplicates are distinct clusters
  
  
  ### RESUME HERE. inefficient.
  
  dat.boot <- bind_rows(lapply(seq_along(samp_sub), function(j) {
    tmp <- by_sub[[samp_sub[j]]]
    tmp$subclass_boot <- paste0(samp_sub[j],".",j)  # new cluster id
    tmp
  }))
  
  # run test neg model
  eqc_itt_fit1 <- coxph(
    Surv(Y3_itt_t_trunc, Y3_itt_factor == "Test Negative") ~ treatment +
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
  eqc_itt_fit2 <- coxph(
    Surv(Y3_itt_t_trunc, Y3_itt_factor == "Test Positive") ~ treatment +
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
    treatHR=exp(eqc_itt_fit2$coefficients["treatment"]) / exp(eqc_itt_fit1$coefficients["treatment"]),
    fluvaxHR=exp(eqc_itt_fit2$coefficients["flu_vax"]) / exp(eqc_itt_fit1$coefficients["flu_vax"])
  ))
  
})

boot.mat <- do.call(rbind, boot.results)






# # Per-Protocol
# 
# 
# eqc_pp_fit1 <- coxph(
#   Surv(Y3_pp_t_trunc, Y3_pp_factor == "Test Negative") ~ treatment +
#     # demog
#     sex_admin + age_years + bmi + race + charlson_cat_fac +
#     # other
#     ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
#     # NEC
#     flu_vax + 
#     cluster(subclass),
#   data = data_Y3
# )
# 
# 
# eqc_pp_fit2 <- coxph(
#   Surv(Y3_pp_t_trunc, Y3_pp_factor == "Test Positive") ~ treatment +
#     # demog
#     sex_admin + age_years + bmi + race + charlson_cat_fac +
#     # other
#     ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
#     # NEC
#     flu_vax + 
#     cluster(subclass),
#   data = data_Y3
# )
# 
# eqc_Y3_cox_pp_tidy <- tibble::tibble(
#   term = c("treatment", "flu_vax"),
#   estimate = c(
#     exp(eqc_pp_fit2$coefficients[1]) / exp(eqc_pp_fit1$coefficients[1]),
#                exp(eqc_pp_fit2$coefficients[6]) / exp(eqc_pp_fit1$coefficients[6])
#     ),
#   std.error = c(NA,NA),
#   statistic = c(NA,NA),
#   p.value = c(NA,NA),
#   conf.low = c(NA,NA),
#   conf.high = c(NA,NA),
# )
# 
# eqc_Y3_cox_pp_tidy
# # write_rds(eqc_Y3_cox_pp_tidy, file = "results/eqc_Y3_cox_pp_tidy.rds") # 2025-12-10


# TND Comparison ----------------------------------------------------------


tnd_fit <- geeglm(
  Y3_itt_factor == "Test Positive" ~ treatment + sex_admin + age_years + ndi +
    bmi + flu_vax + prior_inf + tests_count +
    race + service_region + last_vax_infect_weeks + charlson_cat_fac,
  data = subset(data_Y3, Y3_itt_factor != "Censor"),
  id = subclass,
  family = binomial(link = "logit")
)

tnd_tidy <- tidy(tnd_fit, conf.int = TRUE, exponentiate = TRUE)

tnd_tidy |> print(n=100)
# write_rds(tnd_tidy, file = "results/tnd_tidy.rds") # 2025-12-10


