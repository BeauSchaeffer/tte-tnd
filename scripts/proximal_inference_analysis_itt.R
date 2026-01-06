##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Proximal Inference ITT Analysis

## in progress

# Packages ----------------------------------------------------------------


library(tidyverse)
# library(gt)
# library(gtsummary)
library(tidycmprsk)
library(survival)
# library(ggsurvfit)
# library(riskRegression)


# Data --------------------------------------------------------------------


data_Y3 <- read_rds("cleaned_data/data_Y3.rds")


# Intention to Treat ------------------------------------------------------


# Stage 1

data_Y3 <- data_Y3 |> 
  mutate(Y3_itt_factor = case_when(
    Y3_itt_trunc==0 ~ "Censor",
    Y3_itt_trunc==1 ~ "Test Negative",
    Y3_itt_trunc==2 ~ "Test Positive"
  )) |> 
  mutate(Y3_itt_factor = factor(Y3_itt_factor, levels = c("Censor", "Test Negative", "Test Positive")))

  ### same as ecq_itt_fit1

prox_itt_s1 <- coxph(
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

# Predictions

data_Y3$p <- predict(prox_itt_s1, newdata = data_Y3)

# Stage 2

prox_itt_s2 <- coxph(
  Surv(Y3_itt_t_trunc, Y3_itt_factor == "Test Positive") ~ treatment +
    # demog
    sex_admin + age_years + bmi + race + charlson_cat_fac +
    # other
    ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
    # predictions from S1
    p +
    # NO NEC
    cluster(subclass),
  data = data_Y3
)

prox_itt_tidy <- tidy(prox_itt_s2, conf.int = TRUE, exponentiate = TRUE)
prox_itt_tidy

# write_rds(prox_itt_tidy, file = "results/prox_itt_tidy.rds") # 2025-12-10



# Per-Protocol ------------------------------------------------------------

data_Y3 <- data_Y3 |> 
  mutate(Y3_pp_factor = case_when(
    Y3_pp_trunc==0 ~ "Censor",
    Y3_pp_trunc==1 ~ "Test Negative",
    Y3_pp_trunc==2 ~ "Test Positive"
  )) |> 
  mutate(Y3_pp_factor = factor(Y3_pp_factor, levels = c("Censor", "Test Negative", "Test Positive")))


# Stage 1

prox_pp_s1 <- coxph(
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

# Predictions

data_Y3$p_pp <- predict(prox_pp_s1, newdata = data_Y3)

# Stage 2

prox_pp_s2 <- coxph(
  Surv(Y3_pp_t_trunc, Y3_pp_factor == "Test Positive") ~ treatment + 
    # demog
    sex_admin + age_years + bmi + race + charlson_cat_fac +
    # other
    ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
    # predictions from S1
    p_pp +
    # NO NEC
    cluster(subclass),
  data = data_Y3
)

prox_pp_tidy <- tidy(prox_pp_s2, conf.int = TRUE, exponentiate = TRUE)
# write_rds(prox_pp_tidy, file = "results/prox_pp_tidy.rds")  # 2025-12-10
prox_pp_tidy
prox_itt_tidy

#   -----------------------------------------------------------------------


# -------------------------------------------------------------------------


# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
