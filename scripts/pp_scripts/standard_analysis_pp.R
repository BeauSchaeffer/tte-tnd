##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Standard Analysis
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


data_Y2 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch/data_Y2_weekmatch.rds")

res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results_weekmatch_pp/"


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

summary(std_Y2_cox_pp)

std_Y2_cox_pp_tidy <- tidycmprsk::tidy(std_Y2_cox_pp, conf.int = TRUE, exponentiate = TRUE)
std_Y2_cox_pp_tidy
saveRDS(std_Y2_cox_pp_tidy, paste0(res_path,"std.pp.cox.pointest.rds")) # 2026-06-30









