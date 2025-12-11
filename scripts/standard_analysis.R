##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Standard Analysis

# Packages ----------------------------------------------------------------


library(tidyverse)
library(survival)
library(survminer)
library(tidycmprsk)
library(ggsurvfit)


# Data --------------------------------------------------------------------


data_Y2 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_Y2.rds")


# Kaplan Meier ------------------------------------------------------------


# Intention to Treat

std_Y2_surv_itt <- Surv(time = data_Y2$Y2_itt_t_trunc, event = data_Y2$Y2_itt_trunc)
std_Y2_km_itt <- survfit(std_Y2_surv_itt ~ treatment, data = data_Y2)

std_Y2_km_itt_p <- ggsurvplot(
  std_Y2_km_itt,
  data = data_Y2,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = TRUE,
  xlab = "Weeks",
  ylab = "Survival probability",
  ylim = c(0.90, 1),
  legend.title = "Treatment",
  legend.labs = c("No Booster", "Booster"),
  # palette = "Set1",
  palette = c("#006663", "#FF6B1A"),
  censor = TRUE,
  censor.shape = "+",
  censor.size = 3,
  title = "Kaplan-Meier",
  subtitle = "Intention to Treat (ITT)"
)

std_Y2_km_itt_p

# png("results/std_Y2_km_itt_p.png", width = 8, height = 6, units = "in", res = 300)
# print(std_Y2_km_itt_p)
# dev.off()


# Per-Protocol

std_Y2_surv_pp <- Surv(time = data_Y2$Y2_pp_t_trunc, event = data_Y2$Y2_pp_trunc)
std_Y2_km_pp <- survfit(std_Y2_surv_pp ~ treatment, data = data_Y2)

std_Y2_km_pp_p <- ggsurvplot(
  std_Y2_km_pp,
  data = data_Y2,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = TRUE,
  xlab = "Weeks",
  ylab = "Survival probability",
  ylim = c(0.90, 1),
  legend.title = "Treatment",
  legend.labs = c("No Booster", "Booster"),
  # palette = "Set1",
  palette = c("#006663", "#FF6B1A"),
  censor = TRUE,
  censor.shape = "+",
  censor.size = 3,
  title = "Kaplan-Meier",
  subtitle = "Per-Protocol (PP)"
)

std_Y2_km_pp_p

# png("results/std_Y2_km_pp_p.png", width = 8, height = 6, units = "in", res = 300)
# print(std_Y2_km_pp_p)
# dev.off()


# Cumulative Incidence ----------------------------------------------------


# Intention to Treat

data_Y2 <- data_Y2 |> 
  mutate(Y2_factor = case_when(
    Y2_itt_trunc==0 ~ "Censor",
    Y2_itt_trunc==1 ~ "Test Positive",
  )) |> 
  mutate(Y2_factor = factor(Y2_factor, levels = c("Censor", "Test Positive")))

std_Y2_cif_itt <- cuminc(
  Surv(time = Y2_itt_t_trunc, event = Y2_factor) ~ treatment, 
  data = data_Y2
)

std_Y2_ci_itt_p <- ggcuminc(
  std_Y2_cif_itt, 
  outcome = c("Test Positive"),
  linewidth = 0.5) +
  labs(
    x = "Time",
    y = "Cumulative Incidence",
    title = "Cumulative Incidence",
    subtitle = "Intention to Treat"
  ) +
  add_confidence_interval() +
  scale_color_manual(values = c("0" = "#006663", "1" = "#FF6B1A"),
                     labels = c("0" = "No Booster", "1" = "Booster")) +
  scale_fill_manual(values = c("0" = "#006663", "1" = "#FF6B1A"),
                    labels = c("0" = "No Booster", "1" = "Booster"))

std_Y2_ci_itt_p

# ggsave("results/std_Y2_ci_itt_p.png",std_Y2_ci_itt_p)


# Per-Protocol

std_Y2_cif_pp <- cuminc(
  Surv(time = Y2_pp_t_trunc, event = Y2_factor) ~ treatment, 
  data = data_Y2
)

std_Y2_ci_pp_p <- ggcuminc(
  std_Y2_cif_pp, 
  outcome = c("Test Positive"),
  linewidth = 0.5) +
  labs(
    x = "Time",
    y = "Cumulative Incidence",
    title = "Cumulative Incidence",
    subtitle = "Per-Protocol"
  ) +
  add_confidence_interval() +
  scale_color_manual(values = c("0" = "#006663", "1" = "#FF6B1A"),
                     labels = c("0" = "No Booster", "1" = "Booster")) +
  scale_fill_manual(values = c("0" = "#006663", "1" = "#FF6B1A"),
                    labels = c("0" = "No Booster", "1" = "Booster"))

std_Y2_ci_pp_p

# ggsave("results/std_Y2_ci_pp_p.png",std_Y2_ci_pp_p)


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

  #                               exp(coef) exp(-coef) lower .95 upper .95
  # treatment                     1.1083     0.9023    1.0804    1.1370
  # flu_vax                       1.2855     0.7779    1.2485    1.3235


std_Y2_cox_itt_tidy <- tidycmprsk::tidy(std_Y2_cox_itt, conf.int = TRUE, exponentiate = TRUE)
std_Y2_cox_itt_tidy
write_rds(std_Y2_cox_itt_tidy, file = "results/std_Y2_cox_itt_tidy.rds") # 2025-12-11


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


  #                               exp(coef) exp(-coef) lower .95 upper .95
  # treatment                     0.9664     1.0348    0.9426    0.9908
  # flu_vax                       1.3369     0.7480    1.2990    1.3759

std_Y2_cox_pp_tidy <- tidycmprsk::tidy(std_Y2_cox_pp, conf.int = TRUE, exponentiate = TRUE)
std_Y2_cox_pp_tidy
write_rds(std_Y2_cox_pp_tidy, file = "results/std_Y2_cox_pp_tidy.rds") # 2025-12-11








