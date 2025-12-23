##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Equi Confounding Exploratory

# Packages ----------------------------------------------------------------


library(tidyverse)
library(tidycmprsk)
library(survival)
library(ggsurvfit)
library(riskRegression)
library(geepack)


# Data --------------------------------------------------------------------


data_Y3 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_Y3.rds")


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







