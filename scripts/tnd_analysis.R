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


data_Y3 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch/data_Y3_weekmatch.rds")

data_Y3 <- data_Y3 |> 
  mutate(Y3_itt_factor = case_when(
    Y3_itt_trunc==0 ~ "Censor",
    Y3_itt_trunc==1 ~ "Test Negative",
    Y3_itt_trunc==2 ~ "Test Positive"
  )) |> 
  mutate(Y3_itt_factor = factor(Y3_itt_factor, levels = c("Censor", "Test Negative", "Test Positive")),
         subclass=as.character(subclass))

res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results_weekmatch/"


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

saveRDS(tnd.pointest, paste0(res_path,"tnd.pointest.rds")) # 2026-04-03


# Visualize differences in testing ----------------------------------------

eqc_Y3_cif_itt <- cuminc(
  Surv(Y3_itt_t_trunc, Y3_itt_factor) ~ treatment,
  data = data_Y3
)


# png("figures_draft_wm/tnd.testbehav.plot.png", width = 2400, height=1800, res=300)
par(mar = c(5.1, 5.5, 4.1, 2.1))
plot(NULL,
     xlim = range(c(0, eqc_Y3_cif_itt$tidy$time)),
     ylim = range(c(0, 0.15)),
     xlab="Weeks",
     ylab="Risk",
     main="Nonparametric Risk Curves",
     cex.axis = 1.5,
     cex.lab = 1.5,
     cex.main=1.4
)
mtext("Testing Behavior", side = 3, line = 0.5, font = 3, cex=1.2)
grid()
lines(c(eqc_Y3_cif_itt$tidy$time[eqc_Y3_cif_itt$tidy$outcome=="Test Positive" & eqc_Y3_cif_itt$tidy$strata==0]),
      c(eqc_Y3_cif_itt$tidy$estimate[eqc_Y3_cif_itt$tidy$outcome=="Test Positive" & eqc_Y3_cif_itt$tidy$strata==0]),
      col='#006663', lty=1, lwd=2)
lines(c(eqc_Y3_cif_itt$tidy$time[eqc_Y3_cif_itt$tidy$outcome=="Test Positive" & eqc_Y3_cif_itt$tidy$strata==1]),
      c(eqc_Y3_cif_itt$tidy$estimate[eqc_Y3_cif_itt$tidy$outcome=="Test Positive" & eqc_Y3_cif_itt$tidy$strata==1]),
      col='#FF6B1A', lty=1, lwd=2)
lines(c(eqc_Y3_cif_itt$tidy$time[eqc_Y3_cif_itt$tidy$outcome=="Test Negative" & eqc_Y3_cif_itt$tidy$strata==0]),
      c(eqc_Y3_cif_itt$tidy$estimate[eqc_Y3_cif_itt$tidy$outcome=="Test Negative" & eqc_Y3_cif_itt$tidy$strata==0]),
      col='#006663', lty=2, lwd=2)
lines(c(eqc_Y3_cif_itt$tidy$time[eqc_Y3_cif_itt$tidy$outcome=="Test Negative" & eqc_Y3_cif_itt$tidy$strata==1]),
      c(eqc_Y3_cif_itt$tidy$estimate[eqc_Y3_cif_itt$tidy$outcome=="Test Negative" & eqc_Y3_cif_itt$tidy$strata==1]),
      col='#FF6B1A', lty=2, lwd=2)
legend("topleft",
       legend = c("No Booster", "Booster"),
       col = c('#006663', '#FF6B1A'),
       lty = 1, lwd = 2, cex=1.2,
       bty = "n")
legend("topright",
       legend = c("Test Positive", "Test Negative"),
       col = 'black',
       lty = c(1,2), lwd = 2, cex=1.2,
       bty = "n")
# dev.off()

