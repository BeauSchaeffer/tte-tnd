##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Test Negative Design Analysis ITT
##----- last updated 2026-07-09

# Packages ----------------------------------------------------------------


library(tidyverse)
library(geepack)
library(tidycmprsk)
# library(riskRegression)
library(ggsurvfit)

# Data --------------------------------------------------------------------


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


# TND ---------------------------------------------------------------------


data_Y3_sub <- data_Y3 |>
  filter(Y3_itt_factor != "Censor") |>
  arrange(subclass)

data_Y3_sub <- data_Y3_sub |>
  mutate(subclass_int = as.integer(factor(subclass)))

tnd_fit <- geeglm(
  Y3_itt_factor == "Test Positive" ~ treatment + sex_admin + age_years + ndi +
    bmi + flu_vax + prior_inf + tests_count +
    race + service_region + last_vax_infect_weeks + charlson_cat_fac,
  data = data_Y3_sub,
  id = subclass_int,
  family = binomial(link = "logit")
)

tnd_treatment_coef <- tnd_fit$coefficients["treatment"] |> unname()
tnd_treatment_se <- summary(tnd_fit)$coefficients["treatment", "Std.err"]
tnd_flu_coef <- tnd_fit$coefficients["flu_vax"] |> unname()
tnd_flu_se <- summary(tnd_fit)$coefficients["flu_vax", "Std.err"]


tnd.pointest <- data.frame(
  term = c("treatment", "flu_vax"),
  OR   = c(exp(tnd_treatment_coef), exp(tnd_flu_coef)),
  LCL  = c(exp(tnd_treatment_coef - 1.96*tnd_treatment_se),
           exp(tnd_flu_coef - 1.96*tnd_flu_se)),
  UCL  = c(exp(tnd_treatment_coef + 1.96*tnd_treatment_se),
           exp(tnd_flu_coef + 1.96*tnd_flu_se)))

saveRDS(tnd.pointest, paste0(res_path,"tnd.itt.pointest.rds"))



