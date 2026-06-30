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

res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results_weekmatch_ppp/"


# Cox Model ---------------------------------------------------------------


# Per-Protocol

std_Y2_cox_ppp <- coxph(
  Surv(Y2_pp_t_trunc, Y2_pp_trunc) ~ treatment + 
    sex_admin + age_years + bmi + race + charlson_cat_fac +
    # other
    ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
    # NEC
    flu_vax +
    cluster(subclass),
  data = data_Y2)

summary(std_Y2_cox_ppp)

std_Y2_cox_ppp_tidy <- tidycmprsk::tidy(std_Y2_cox_ppp, conf.int = TRUE, exponentiate = TRUE)
std_Y2_cox_ppp_tidy
saveRDS(std_Y2_cox_ppp_tidy, paste0(res_path,"std.ppp.cox.pointest.rds")) # 2026-06-30









