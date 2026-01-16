##----- Kaiser Causal TTE-TND
##----- Combine Results


# Packages ----------------------------------------------------------------


library(tidyverse)


# Data --------------------------------------------------------------------


res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results/"

# STD Cox

  std.itt.cox.pointest <- readRDS(paste0(res_path,"std.itt.cox.pointest.rds"))

# STD Pooled

  # std_pooled_itt_model <- readRDS(paste0(res_path, "std_pooled_itt_model.rds"))
  std.itt.risk.pointest <- readRDS(paste0(res_path, "std.itt.risk.pointest.rds"))
  std.itt.boot.long <- readRDS(paste0(res_path, "std.itt.boot.long.rds"))

# TND

  tnd.pointest <- readRDS(paste0(res_path,"tnd.pointest.rds"))

# EQC Cox

  eqc.itt.cox.pointest <- readRDS(paste0(res_path, "eqc.itt.cox.pointest.rds"))
  eqc.itt.cox.boot.long <- readRDS(paste0(res_path, "eqc.itt.cox.boot.long.rds"))

# EQC Pooled
  
  eqc.itt.risk.pointest <- readRDS(paste0(res_path, "eqc.itt.risk.pointest.rds"))
  eqc.itt.boot.long <- readRDS(paste0(res_path, "eqc.itt.boot.long.rds"))
  
# PCI Cox
  
  pci.itt.cox.pointest <- readRDS(paste0(res_path, "pci.itt.cox.pointest.rds"))
  pci.itt.cox.boot.long <- readRDS(paste0(res_path, "pci.itt.cox.boot.long.rds"))
  
# PCI Pooled
  pci.itt.risk.pointest <- readRDS(paste0(res_path, "pci.itt.risk.pointest.rds"))
  
  pci_rep_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results/pci_boot_reps"
  
  pci_rep_files <- list.files(
    pci_rep_path,
    pattern = "^pci_itt_boot_rep_\\d{3}\\.rds$",
    full.names = TRUE
  )
  
  pci.itt.boot.long <- pci_rep_files |>
    lapply(readRDS) |>
    lapply(\(m) as.data.frame(m)) |>
    bind_rows() |>
    as_tibble() |>
    mutate(
      sim = as.integer(sim),
      time_end = as.integer(time_end),
      risk0 = as.numeric(risk0),
      risk1 = as.numeric(risk1)
    ) |>
    arrange(sim, time_end)
  
  saveRDS(pci.itt.boot.long, file.path(res_path, "pci.itt.boot.long.rds")) # 2025-01-16
  
  
  
# Boot CI functions -------------------------------------------------------

  
pooled.boot.ci <- function(point.est, boot.long, alpha = 0.05){
    
    # decide which column to treat as "risk0"
    risk0_var <- if ("risk0corr" %in% names(boot.long)) "risk0corr" else "risk0"
    
    boot.ci <- boot.long |>
      group_by(time_end) |>
      summarise(
        risk0_lo = quantile(.data[[risk0_var]], probs = alpha/2, na.rm = TRUE),
        risk0_hi = quantile(.data[[risk0_var]], probs = 1 - alpha/2, na.rm = TRUE),
        risk1_lo = quantile(risk1, probs = alpha/2, na.rm = TRUE),
        risk1_hi = quantile(risk1, probs = 1 - alpha/2, na.rm = TRUE),
        .groups = "drop"
      )
    
    # ensure point estimate uses the same risk0 definition
    point.est.use <- point.est |>
      mutate(risk0 = if (risk0_var == "risk0corr") risk0corr else risk0)
    
    boot.point.ci <- boot.ci |>
      left_join(point.est.use, by = "time_end") |>
      select(-sim, -any_of("risk0corr")) |>
      relocate(time_end, risk0, risk0_lo, risk0_hi, risk1, risk1_lo, risk1_hi)
    
    return(boot.point.ci)
}

eqc.cox.boot.ci <- function(point.est, boot.long, alpha = 0.05){

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

pci.cox.boot.ci <- function(point.est, boot.long, alpha = 0.05){
  
  boot.ci <- boot.long |>
    summarise(
      treatHR_lo = quantile(treatHR, probs = alpha/2, na.rm = TRUE),
      treatHR_hi = quantile(treatHR, probs = 1 - alpha/2, na.rm = TRUE)
    )
  
  treatHR <- point.est |>
    filter(term == "treatment") |>
    pull(estimate) |>
    first()
  
  point.df <- data.frame(treatHR = treatHR)
  
  
  bind_cols(point.df, boot.ci) |>
    relocate(treatHR, treatHR_lo, treatHR_hi)
  
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
# saveRDS(std.itt.risks.ci, paste0(res_path, "std.itt.risks.ci.rds")) # 2025-12-23

# png("figures_draft/std.itt.risks.ci.plot.png", width = 2400, height=1800, res=300)
plot.risk.with.boot.ci(std.itt.risks.ci, title.sub  = "Standard TTE (ITT)")
# dev.off()


# EQC ITT -----------------------------------------------------------------


eqc.itt.HRs.ci <- eqc.cox.boot.ci(eqc.itt.cox.pointest, eqc.itt.cox.boot.long)
# saveRDS(eqc.itt.HRs.ci, paste0(res_path, "eqc.itt.HRs.ci.rds")) # 2025-12-30

eqc.itt.risks.ci <- pooled.boot.ci(point.est = eqc.itt.risk.pointest, boot.long = eqc.itt.boot.long)
# saveRDS(eqc.itt.risks.ci, paste0(res_path, "eqc.itt.risks.ci.rds")) # 2026-01-06

# png("figures_draft/eqc.itt.risks.ci.plot.png", width = 2400, height=1800, res=300)
plot.risk.with.boot.ci(eqc.itt.risks.ci, title.sub  = "Equi-confounding (ITT)")
# dev.off()


# PCI ITT -----------------------------------------------------------------

pci.itt.HRs.ci <- pci.cox.boot.ci(pci.itt.cox.pointest, pci.itt.cox.boot.long)
# saveRDS(pci.itt.HRs.ci, paste0(res_path, "pci.itt.HRs.ci.rds")) # 2026-01-12

pci.itt.risks.ci <- pooled.boot.ci(point.est = pci.itt.risk.pointest, boot.long = pci.itt.boot.long)
# saveRDS(pci.itt.risks.ci, paste0(res_path, "pci.itt.risks.ci.rds")) # 2026-01-16

# png("figures_draft/pci.itt.risks.ci.plot.png", width = 2400, height=1800, res=300)
plot.risk.with.boot.ci(pci.itt.risks.ci, title.sub  = "Proximal inference (ITT)")
# dev.off()





