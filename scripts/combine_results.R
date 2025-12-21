##----- Kaiser Causal TTE-TND
##----- Combine Results

# Data --------------------------------------------------------------------

res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results/"

std.itt.boot.long <- readRDS(paste0(res_path, "std.itt.boot.long.rds"))

std.itt.risk.pointest <- readRDS(paste0(res_path, "std.itt.risk.pointest.rds"))


# STD ITT -----------------------------------------------------------------


# Pointwise bootstrap CIs (95%)
alpha <- 0.05
std.itt.boot.ci <- std.itt.boot.long |>
  dplyr::group_by(time_end) |>
  dplyr::summarise(
    risk0_lo = stats::quantile(risk0, probs = alpha/2, na.rm = TRUE),
    risk0_hi = stats::quantile(risk0, probs = 1 - alpha/2, na.rm = TRUE),
    risk1_lo = stats::quantile(risk1, probs = alpha/2, na.rm = TRUE),
    risk1_hi = stats::quantile(risk1, probs = 1 - alpha/2, na.rm = TRUE),
    .groups = "drop"
  )

# Plot
par(mar = c(5.1, 5.5, 4.1, 2.1))
plot(NULL,
     xlim = range(c(0, std.itt.risk.pointest$time_end, std.itt.boot.ci$time_end), na.rm = TRUE),
     ylim = c(0, 0.075),
     xlab="Weeks",
     ylab="Risk",
     main="Risk Curves",
     cex.axis = 1.5,
     cex.lab = 1.5,
     cex.main=1.4
)
mtext("Standard TTE (ITT)", side = 3, line = 0.5, font = 3, cex=1.2)
grid()

# Bootstrap CI ribbons (behind lines)
polygon(
  x = c(std.itt.boot.ci$time_end, rev(std.itt.boot.ci$time_end)),
  y = c(std.itt.boot.ci$risk0_lo, rev(std.itt.boot.ci$risk0_hi)),
  # col = adjustcolor("#006663", alpha.f = 0.80),
  border = NA
)
polygon(
  x = c(std.itt.boot.ci$time_end, rev(std.itt.boot.ci$time_end)),
  y = c(std.itt.boot.ci$risk1_lo, rev(std.itt.boot.ci$risk1_hi)),
  # col = adjustcolor("#FF6B1A", alpha.f = 0.80),
  border = NA
)

# Point estimate lines (assumes point estimate has columns time_end, risk0, risk1)
lines(c(0, std.itt.risk.pointest$time_end),
      c(0, std.itt.risk.pointest$risk0),
      col = "#006663", lty = 1, lwd = 1)

lines(c(0, std.itt.risk.pointest$time_end),
      c(0, std.itt.risk.pointest$risk1),
      col = "#FF6B1A", lty = 1, lwd = 1)

legend("topleft",
       legend = c("No Booster", "Booster"),
       col = c("#006663", "#FF6B1A"),
       lty = 1, lwd = 4, cex=1.2,
       bty = "n")

