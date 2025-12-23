##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Equi Confounding Analysis

# setwd("~/Desktop/Research/Kaiser/KP_analysis")


# Packages ----------------------------------------------------------------


library(tidyverse)
# library(gt)
# library(gtsummary)
library(tidycmprsk)
library(survival)
library(ggsurvfit)
library(riskRegression)
library(geepack)


# Data --------------------------------------------------------------------


data_Y3 <- readr::read_rds("cleaned_data/data_Y3.rds")


# Cumulative Incidence ----------------------------------------------------


# Intention to Treat

data_Y3 <- data_Y3 |> 
  mutate(Y3_itt_factor = case_when(
    Y3_itt_trunc==0 ~ "Censor",
    Y3_itt_trunc==1 ~ "Test Negative",
    Y3_itt_trunc==2 ~ "Test Positive"
  )) |> 
  mutate(Y3_itt_factor = factor(Y3_itt_factor, levels = c("Censor", "Test Negative", "Test Positive")))

eqc_Y3_cif_itt <- cuminc(
  Surv(Y3_itt_t_trunc, Y3_itt_factor) ~ treatment, 
  data = data_Y3
)

eqc_Y3_ci_itt_p <- ggcuminc(
  eqc_Y3_cif_itt,
  aes(color = group, linetype = outcome),
  outcome = c("Test Negative", "Test Positive"),
  linewidth = 0.5) +
  labs(
    x = "Time",
    y = "Cumulative incidence",
    title = "Cumulative incidence",
    subtitle = "Intention to Treat",
    color = "Treatment",
    linetype = "Outcome"
  ) +
  add_confidence_interval() +
  scale_color_manual(
    values = c("0" = "#006663", "1" = "#FF6B1A"),
    labels = c("0" = "No Booster", "1" = "Booster")
  ) +
  scale_fill_manual(
    values = c("0" = "#006663", "1" = "#FF6B1A"),
    guide = "none"
  ) +
  scale_linetype_manual(
    values = c("Test Negative" = "dashed", "Test Positive" = "solid")
  ) +
  ylim(0, 0.20)

eqc_Y3_ci_itt_p

# ggsave("results/eqc_Y3_ci_itt_p.png",eqc_Y3_ci_itt_p) # 2025-12-10

# Per-Protocol

data_Y3 <- data_Y3 |> 
  mutate(Y3_pp_factor = case_when(
    Y3_pp_trunc==0 ~ "Censor",
    Y3_pp_trunc==1 ~ "Test Negative",
    Y3_pp_trunc==2 ~ "Test Positive"
  )) |> 
  mutate(Y3_pp_factor = factor(Y3_pp_factor, levels = c("Censor", "Test Negative", "Test Positive")))

eqc_Y3_cif_pp <- cuminc(
  Surv(Y3_pp_t_trunc, Y3_pp_factor) ~ treatment, 
  data = data_Y3
)

eqc_Y3_ci_pp_p <- ggcuminc(
  eqc_Y3_cif_pp,
  aes(color = group, linetype = outcome),
  outcome = c("Test Negative", "Test Positive"),
  linewidth = 0.5) +
  labs(
    x = "Time",
    y = "Cumulative incidence",
    title = "Cumulative incidence",
    subtitle = "Per-Protocol",
    color = "Treatment",
    linetype = "Outcome"
  ) +
  add_confidence_interval() +
  scale_color_manual(
    values = c("0" = "#006663", "1" = "#FF6B1A"),
    labels = c("0" = "No Booster", "1" = "Booster")
  ) +
  scale_fill_manual(
    values = c("0" = "#006663", "1" = "#FF6B1A"),
    guide = "none"
  ) +
  scale_linetype_manual(
    values = c("Test Negative" = "dashed", "Test Positive" = "solid")
  ) +
  ylim(0, 0.20)

eqc_Y3_ci_pp_p

# ggsave("results/eqc_Y3_ci_pp_p.png",eqc_Y3_ci_pp_p) # 2025-12-10


# EQ approach -------------------------------------------------------------


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

eqc_Y3_cox_itt_tidy <- tibble::tibble(
  term = c("treatment", "flu_vax"),
  estimate = c(
    exp(eqc_itt_fit2$coefficients[1]) / exp(eqc_itt_fit1$coefficients[1]),
               exp(eqc_itt_fit2$coefficients[6]) / exp(eqc_itt_fit1$coefficients[6])
    ),
  std.error = c(NA,NA),
  statistic = c(NA,NA),
  p.value = c(NA,NA),
  conf.low = c(NA,NA),
  conf.high = c(NA,NA),
)

eqc_Y3_cox_itt_tidy
# write_rds(eqc_Y3_cox_itt_tidy, file = "results/eqc_Y3_cox_itt_tidy.rds") # 2025-12-10


# Per-Protocol


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

eqc_Y3_cox_pp_tidy <- tibble::tibble(
  term = c("treatment", "flu_vax"),
  estimate = c(
    exp(eqc_pp_fit2$coefficients[1]) / exp(eqc_pp_fit1$coefficients[1]),
               exp(eqc_pp_fit2$coefficients[6]) / exp(eqc_pp_fit1$coefficients[6])
    ),
  std.error = c(NA,NA),
  statistic = c(NA,NA),
  p.value = c(NA,NA),
  conf.low = c(NA,NA),
  conf.high = c(NA,NA),
)

eqc_Y3_cox_pp_tidy
# write_rds(eqc_Y3_cox_pp_tidy, file = "results/eqc_Y3_cox_pp_tidy.rds") # 2025-12-10


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


