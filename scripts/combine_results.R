##----- Kaiser Causal TTE-TND
##----- Combine Results


# Packages ----------------------------------------------------------------


library(tidyverse)


# Data --------------------------------------------------------------------


res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results/"

# STD Cox - complete

  std.itt.cox.pointest <- readRDS(paste0(res_path,"std.itt.cox.pointest.rds"))

# STD Pooled - rerunning with paired sampling

  # std_pooled_itt_model <- readRDS(paste0(res_path, "std_pooled_itt_model.rds"))
  std.itt.risk.pointest <- readRDS(paste0(res_path, "std.itt.risk.pointest.rds"))
  std.itt.boot.long <- readRDS(paste0(res_path, "std.itt.boot.long.rds"))

# TND - complete

  tnd.pointest <- readRDS(paste0(res_path,"tnd.pointest.rds"))

# EQC Cox - complete

  eqc.itt.cox.pointest <- readRDS(paste0(res_path, "eqc.itt.cox.pointest.rds"))
  eqc.itt.cox.boot.long <- readRDS(paste0(res_path, "eqc.itt.cox.boot.long.rds"))

# EQC Pooled - incomplete
  
# PCI Cox - incomplete
  
# PCI Pooled - incomplete
  
  
# Boot CI functions -------------------------------------------------------

pooled.boot.ci <- function(point.est, boot.long, alpha = 0.05){
  
  boot.ci <- boot.long |>
    group_by(time_end) |>
    summarise(
      risk0_lo = quantile(risk0, probs = alpha/2, na.rm = TRUE),
      risk0_hi = quantile(risk0, probs = 1 - alpha/2, na.rm = TRUE),
      risk1_lo = quantile(risk1, probs = alpha/2, na.rm = TRUE),
      risk1_hi = quantile(risk1, probs = 1 - alpha/2, na.rm = TRUE),
      .groups = "drop"
    )
  
  boot.point.ci <- boot.ci |>
    left_join(point.est, by = "time_end") |>
    select(-sim) |>
    relocate(time_end, risk0, risk0_lo, risk0_hi, risk1, risk1_lo, risk1_hi)
  
  return(boot.point.ci)
}
  
cox.boot.ci <- function(point.est, boot.long, alpha = 0.05){
  
  boot.ci <- boot.long |> 
    summarise(
      treatHR_lo = quantile(treatHR, probs = alpha/2, na.rm = TRUE),
      treatHR_hi = quantile(treatHR, probs = 1 - alpha/2, na.rm = TRUE),
      fluvaxHR_lo = quantile(fluvaxHR, probs = alpha/2, na.rm = TRUE),
      fluvaxHR_hi = quantile(fluvaxHR, probs = 1 - alpha/2, na.rm = TRUE)
    )
  
  point.df <- data.frame(
    treatHR = unname(eqc.itt.cox.pointest["treatHR"]),
    fluvaxHR = unname(eqc.itt.cox.pointest["fluvaxHR"])
  )
  
  bind_cols(point.df, boot.ci) |> 
    relocate(treatHR, treatHR_lo, treatHR_hi, fluvaxHR, fluvaxHR_lo, fluvaxHR_hi)
  
}

# Plotting function -------------------------------------------------------

plot.risk.with.boot.ci <- function(risks.and.cis,
                                   title.main = "Risk Curves",
                                   title.sub  = "Approach (Effect)",
                                   xlab = "Weeks",
                                   ylab = "Risk",
                                   ylim = c(0, 0.075),
                                   col0 = "#006663",
                                   col1 = "#FF6B1A",
                                   ribbon.alpha = 0.25,
                                   lwd.lines = 1,
                                   lwd.legend = 4,
                                   legend_pos = "topleft",
                                   add.grid = TRUE,
                                   mar = c(5.1, 5.5, 4.1, 2.1),
                                   cex.axis = 1.5,
                                   cex.lab  = 1.5,
                                   cex.main = 1.4,
                                   cex.sub  = 1.2,
                                   legend.cex = 1.2) {
  
  req_cols <- c("time_end",
                "risk0", "risk0_lo", "risk0_hi",
                "risk1", "risk1_lo", "risk1_hi")
  
  missing <- setdiff(req_cols, names(risks.and.cis))
  if (length(missing) > 0) {
    stop("`risks.and.cis` is missing columns: ", paste(missing, collapse = ", "))
  }
  
  risks.and.cis <- risks.and.cis |> arrange(time_end)

  xlim <- range(c(0, risks.and.cis$time_end), na.rm = TRUE)
  
  par(mar = mar)
  plot(NULL,
       xlim = xlim,
       ylim = ylim,
       xlab = xlab,
       ylab = ylab,
       main = title.main,
       cex.axis = cex.axis,
       cex.lab  = cex.lab,
       cex.main = cex.main)
  
  if (!is.null(title.sub) && nzchar(title.sub)) {
    mtext(title.sub, side = 3, line = 0.5, font = 3, cex = cex.sub)
  }
  if (isTRUE(add.grid)) grid()
  
  # CI ribbons
  polygon(
    x = c(risks.and.cis$time_end, rev(risks.and.cis$time_end)),
    y = c(risks.and.cis$risk0_lo, rev(risks.and.cis$risk0_hi)),
    col = adjustcolor(col0, alpha.f = ribbon.alpha),
    border = NA
  )
  
  polygon(
    x = c(risks.and.cis$time_end, rev(risks.and.cis$time_end)),
    y = c(risks.and.cis$risk1_lo, rev(risks.and.cis$risk1_hi)),
    col = adjustcolor(col1, alpha.f = ribbon.alpha),
    border = NA
  )
  
  # point estimate lines
  lines(c(0, risks.and.cis$time_end),
        c(0, risks.and.cis$risk0),
        col = col0, lty = 1, lwd = lwd.lines)
  
  lines(c(0, risks.and.cis$time_end),
        c(0, risks.and.cis$risk1),
        col = col1, lty = 1, lwd = lwd.lines)

  legend(legend_pos,
         legend = c("No Booster", "Booster"),
         col = c(col0, col1),
         lty = 1, lwd = lwd.legend, cex = legend.cex,
         bty = "n")
  
  invisible(risks.and.cis)
}


# STD ITT -----------------------------------------------------------------

std.itt.risks.ci <- pooled.boot.ci(point.est = std.itt.risk.pointest, boot.long = std.itt.boot.long)
# saveRDS(std_itt_risks_ci, paste0(res_path, "std.itt.risks.ci.rds")) # 2025-12-23

plot.risk.with.boot.ci(std_itt_risks_ci, title.sub  = "Standard TTE (ITT)")


# EQC ITT -----------------------------------------------------------------


eqc.itt.HRs.ci <- cox.boot.ci(eqc.itt.cox.pointest, eqc.itt.cox.boot.long)
# saveRDS(eqc.itt.HRs.ci, paste0(res_path, "eqc.itt.HRs.ci.rds")) # 2025-12-30







