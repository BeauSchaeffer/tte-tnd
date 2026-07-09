##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Standard Analysis
##----- last updated 2026-07-09

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

data_Y3 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch/data_Y3_weekmatch.rds")

data_Y3 <- data_Y3 |> 
  mutate(Y3_itt_factor = case_when(
    Y3_itt_trunc==0 ~ "Censor",
    Y3_itt_trunc==1 ~ "Test Negative",
    Y3_itt_trunc==2 ~ "Test Positive"
  )) |> 
  mutate(Y3_itt_factor = factor(Y3_itt_factor, levels = c("Censor", "Test Negative", "Test Positive")),
         subclass=as.character(subclass))

res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results_itt.2/"


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

std_Y2_cox_itt_tidy <- tidycmprsk::tidy(std_Y2_cox_itt, conf.int = TRUE, exponentiate = TRUE)
saveRDS(std_Y2_cox_itt_tidy, paste0(res_path,"std.itt.cox.pointest.rds"))


# HR for NCO

std_cox_nco_itt <- coxph(
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

std_cox_nco_itt_tidy <- tidycmprsk::tidy(std_cox_nco_itt, conf.int = TRUE, exponentiate = TRUE)
saveRDS(std_cox_nco_itt_tidy, paste0(res_path,"std.cox.nco.itt.tidy.rds"))

