##----- Beau Schaeffer
##----- Kaiser TTE-TND
##----- Censoring by arm


# Packages ----------------------------------------------------------------


library(tidyverse)
library(tidycmprsk)


# Data --------------------------------------------------------------------


data_Y3 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch.3/data_Y3_weekmatch.rds")

plot_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/plots_pp.5/"

max_follow <- 52


# Build censoring-as-event data -------------------------------------------

# status:  1 = treatment deviation censoring (event of interest)
#          2 = testing event (competing)
#          0 = administrative censoring / completed follow-up

cens_dat <- as.data.frame(data_Y3) |>
  mutate(
    pp_cens = Y3_pp_t < Y3_itt_t & Y3_pp_t < max_follow,
    
    status = case_when(
      pp_cens                    ~ 1L,
      Y3_pp_trunc %in% c(1L, 2L) ~ 2L,
      TRUE                       ~ 0L
    ),
    status = factor(status, levels = c(0, 1, 2)),
    arm    = factor(treatment, levels = c(0, 1))
  )

# check
# cens_dat |> count(treatment, status)


# Plot --------------------------------------------------------------------


cens.cuminc <- cuminc(Surv(Y3_pp_t_trunc, status) ~ arm, data = cens_dat)

d  <- tidy(cens.cuminc) |> filter(outcome == "1")
d0 <- d |> filter(strata == "0")
d1 <- d |> filter(strata == "1")


png(paste0(plot_path, "censoring.cuminc.pp.png"), width = 2400, height = 1800, res = 300)

par(mar = c(5, 5.5, 4, 2))

plot(NA, xlim = c(0, max_follow), ylim = c(0, 0.5),
     xlab = "Weeks", ylab = "Cumulative incidence of censoring",
     main = "Treatment Deviation Censoring",
     font.main = 1, cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.3, las = 1)
grid(col = "grey85", lty = 1, lwd = 0.5)
polygon(c(d0$time, rev(d0$time)), c(d0$conf.low, rev(d0$conf.high)),
        col = adjustcolor("#006663", alpha.f = 0.2), border = NA)
polygon(c(d1$time, rev(d1$time)), c(d1$conf.low, rev(d1$conf.high)),
        col = adjustcolor("#FF6B1A", alpha.f = 0.2), border = NA)

points(d0$time, d0$estimate, pch = 16, col = "#006663", cex = 0.6)
points(d1$time, d1$estimate, pch = 16, col = "#FF6B1A", cex = 0.6)

legend("topleft", legend = c("No booster", "Booster"),
       col = c("#006663", "#FF6B1A"), pch = 16, bty = "n", cex = 1.1)

dev.off()