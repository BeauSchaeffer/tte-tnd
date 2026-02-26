##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Test Negative Design Analysis

# Packages ----------------------------------------------------------------


library(tidyverse)
library(geepack)
library(tidycmprsk)
# library(riskRegression)
library(ggsurvfit)

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


tnd.pointest <- data.frame(
  term = c("treatment", "flu_vax"),
  OR   = c(exp(tnd_treatment_coef), exp(tnd_flu_coef)),
  LCL  = c(exp(tnd_treatment_coef - 1.96*tnd_treatment_se),
           exp(tnd_flu_coef - 1.96*tnd_flu_se)),
  UCL  = c(exp(tnd_treatment_coef + 1.96*tnd_treatment_se),
           exp(tnd_flu_coef + 1.96*tnd_flu_se)))

saveRDS(tnd.pointest, paste0(res_path,"tnd.pointest.rds")) # 2025-12-26


# Visualize differences in testing ----------------------------------------

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


### base R

par(mar = c(5.1, 5.5, 4.1, 2.1))
plot(NULL,
     xlim = range(c(0, eqc_Y3_cif_itt$tidy$time)),
     ylim = range(c(0, 0.1)),
     xlab="Weeks",
     ylab="Risk",
     main="Nonparametric Risk Curves",
     cex.axis = 1.5,
     cex.lab = 1.5,
     cex.main=1.4
)
# mtext("Standard TTE (ITT)", side = 3, line = 0.5, font = 3, cex=1.2)
grid()
lines(c(0, eqc_Y3_cif_itt$tidy$time[eqc_Y3_cif_itt$tidy$outcome=="Test Positive" & eqc_Y3_cif_itt$tidy$strata==0]), 
      c(0, eqc_Y3_cif_itt$tidy$estimate[eqc_Y3_cif_itt$tidy$outcome=="Test Positive" & eqc_Y3_cif_itt$tidy$strata==0]),
      col='#006663', lty=1, lwd=4)
# lines(c(0, std_itt_A1.long$time_end), c(0, std_itt_A1.long$risk), col='#FF6B1A', lty=1, lwd=4)
# legend("topleft",
#        legend = c("No Booster", "Booster"),
#        col = c('#006663', '#FF6B1A'),
#        lty = 1, lwd = 4, cex=1.2,
#        bty = "n")

