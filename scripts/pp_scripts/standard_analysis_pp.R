##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Standard Analysis
##----- Per-protocol, no censoring weights
##----- last updated 2026-07-16

# Packages ----------------------------------------------------------------


library(tidyverse)
library(tidycmprsk)
library(survival)
library(ggsurvfit)
library(riskRegression)
library(geepack)
library(data.table)

# Data --------------------------------------------------------------------


data_Y2 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch/data_Y2_weekmatch.rds")

data_Y3 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch/data_Y3_weekmatch.rds")

data_Y3 <- data_Y3 |> 
  mutate(Y3_pp_factor = case_when(
    Y3_pp_trunc==0 ~ "Censor",
    Y3_pp_trunc==1 ~ "Test Negative",
    Y3_pp_trunc==2 ~ "Test Positive"
  )) |> 
  mutate(Y3_pp_factor = factor(Y3_pp_factor, levels = c("Censor", "Test Negative", "Test Positive")),
         subclass=as.character(subclass))

res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results_pp.2/"


# Cox Model ---------------------------------------------------------------


# Per-Protocol

std_Y2_cox_pp <- coxph(
  Surv(Y2_pp_t_trunc, Y2_pp_trunc) ~ treatment + 
    sex_admin + age_years + bmi + race + charlson_cat_fac +
    # other
    ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
    # NEC
    flu_vax +
    cluster(subclass),
  data = data_Y2)

std_Y2_cox_pp_tidy <- tidycmprsk::tidy(std_Y2_cox_pp, conf.int = TRUE, exponentiate = TRUE)
saveRDS(std_Y2_cox_pp_tidy, paste0(res_path,"std.pp.cox.pointest.rds"))


# HR for NCO

std_cox_nco_pp <- coxph(
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

std_cox_nco_pp_tidy <- tidycmprsk::tidy(std_cox_nco_pp, conf.int = TRUE, exponentiate = TRUE)
saveRDS(std_cox_nco_pp_tidy, paste0(res_path,"std.cox.nco.pp.tidy.rds"))

