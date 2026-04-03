##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Standard Analysis

# Packages ----------------------------------------------------------------


library(tidyverse)
library(tidycmprsk)
library(survival)
library(ggsurvfit)
library(riskRegression)
library(geepack)
library(data.table)

# Data --------------------------------------------------------------------


# data_Y2 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_Y2.rds")
data_Y2 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch/data_Y2_weekmatch.rds")

res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results_weekmatch/"


# Cox Model ---------------------------------------------------------------


# Intention to Treat

std_Y2_cox_itt <- coxph(
  Surv(Y2_itt_t_trunc, Y2_itt_trunc) ~ treatment + 
    sex_admin + age_years + bmi + race + charlson_cat_fac +
    # other
    ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
    # NEC
    flu_vax +
    cluster(subclass),
  data = data_Y2)

summary(std_Y2_cox_itt)

std_Y2_cox_itt_tidy <- tidycmprsk::tidy(std_Y2_cox_itt, conf.int = TRUE, exponentiate = TRUE)
std_Y2_cox_itt_tidy
saveRDS(std_Y2_cox_itt_tidy, paste0(res_path,"std.itt.cox.pointest.rds")) # 2025-12-26


# # Per-Protocol
# 
# std_Y2_cox_pp <- coxph(
#   Surv(Y2_pp_t_trunc, Y2_pp_trunc) ~ treatment + 
#     sex_admin + age_years + bmi + race + charlson_cat_fac +
#     # other
#     ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
#     # NEC
#     flu_vax +
#     cluster(subclass),
#   data = data_Y2)
# 
# summary(std_Y2_cox_pp)
# 
# 
#   #                               exp(coef) exp(-coef) lower .95 upper .95
#   # treatment                     0.9664     1.0348    0.9426    0.9908
#   # flu_vax                       1.3369     0.7480    1.2990    1.3759
# 
# std_Y2_cox_pp_tidy <- tidycmprsk::tidy(std_Y2_cox_pp, conf.int = TRUE, exponentiate = TRUE)
# std_Y2_cox_pp_tidy
# write_rds(std_Y2_cox_pp_tidy, file = "results/std_Y2_cox_pp_tidy.rds") # 2025-12-11








