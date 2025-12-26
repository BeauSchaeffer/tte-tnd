# ##----- Beau Schaeffer
# ##----- Kaiser Causal TTE-TND
# ##----- Equi Confounding Analysis PP
# 
# # Packages ----------------------------------------------------------------
# 
# 
# library(tidyverse)
# library(tidycmprsk)
# library(survival)
# library(ggsurvfit)
# library(riskRegression)
# library(geepack)
# library(data.table)
# 
# 
# # Data --------------------------------------------------------------------
# 
# 
# data_Y3 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_Y3.rds")
# 
# data_Y3 <- data_Y3 |> 
#   mutate(Y3_itt_factor = case_when(
#     Y3_itt_trunc==0 ~ "Censor",
#     Y3_itt_trunc==1 ~ "Test Negative",
#     Y3_itt_trunc==2 ~ "Test Positive"
#   )) |> 
#   mutate(Y3_itt_factor = factor(Y3_itt_factor, levels = c("Censor", "Test Negative", "Test Positive")),
#          subclass=as.character(subclass))
# 
# res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results/"
# 
# 

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
