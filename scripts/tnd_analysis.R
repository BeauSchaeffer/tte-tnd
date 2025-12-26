##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Test Negative Design Analysis

# Packages ----------------------------------------------------------------


library(tidyverse)
library(geepack)


# Data --------------------------------------------------------------------


data_Y3 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_Y3.rds")

data_Y3 <- data_Y3 |> 
  mutate(Y3_itt_factor = case_when(
    Y3_itt_trunc==0 ~ "Censor",
    Y3_itt_trunc==1 ~ "Test Negative",
    Y3_itt_trunc==2 ~ "Test Positive"
  )) |> 
  mutate(Y3_itt_factor = factor(Y3_itt_factor, levels = c("Censor", "Test Negative", "Test Positive")),
         subclass=as.character(subclass))

res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results/"


# TND ---------------------------------------------------------------------


tnd_fit <- geeglm(
  Y3_itt_factor == "Test Positive" ~ treatment + sex_admin + age_years + ndi +
    bmi + flu_vax + prior_inf + tests_count +
    race + service_region + last_vax_infect_weeks + charlson_cat_fac,
  data = subset(data_Y3, Y3_itt_factor != "Censor"),
  id = subclass,
  family = binomial(link = "logit")
)

tnd_treatment_coef <- tnd_fit$coefficients["treatment"] |> unname()
tnd_treatment_se <- summary(tnd_fit)$coefficients["treatment", "Std.err"]
tnd_flu_coef <- tnd_fit$coefficients["flu_vax"] |> unname()
tnd_flu_se <- summary(tnd_fit)$coefficients["flu_vax", "Std.err"]


tnd.pointest <- cbind(
  treatment=exp(tnd_treatment_coef),
  treatment_lo=exp(tnd_treatment_coef - 1.96*tnd_treatment_se),
  treatment_hi=exp(tnd_treatment_coef + 1.96*tnd_treatment_se),
  
  flu_vax=exp(tnd_flu_coef),
  flu_vax_lo=exp(tnd_flu_coef - 1.96*tnd_flu_se),
  flu_vax_hi=exp(tnd_flu_coef + 1.96*tnd_flu_se),
)

saveRDS(eqc.itt.cox.pointest, paste0(res_path,"eqc.itt.cox.pointest.rds")) # 2025-12-26

tnd_tidy <- tidy(tnd_fit, conf.int = TRUE, exponentiate = TRUE)

tnd_tidy |> print(n=100)
# write_rds(tnd_tidy, file = "results/tnd_tidy.rds") # 2025-12-10