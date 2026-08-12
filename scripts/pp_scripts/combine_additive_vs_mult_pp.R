##----- Kaiser Causal TTE-TND
##----- Additive vs multiplicative marginal CIF -- combine + overlay plot
##----- Per-protocol
##----- last updated 2026-08-07
##-----
##----- Computes bootstrap CIs for the ADDITIVE marginal CIF curves and overlays
##----- them on the existing MULTIPLICATIVE (pooled-logistic) curves for EQC and PCI.


# Packages ----------------------------------------------------------------


library(tidyverse)


# Paths -------------------------------------------------------------------


res_path  <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results_pp.5/"
plot_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/plots_pp.5/"


# Bootstrap CI helper (additive) ------------------------------------------


add.boot.ci <- function(point.est, boot.long, alpha = 0.05){
  # use the corrected curve as risk0 when present (EQC), else the naive risk0 (PCI)
  risk0_var <- if ("risk0corr" %in% names(boot.long)) "risk0corr" else "risk0"
  ci <- boot.long |>
    group_by(time_end) |>
    summarise(
      risk0_lo = quantile(.data[[risk0_var]], alpha/2,     na.rm = TRUE),
      risk0_hi = quantile(.data[[risk0_var]], 1 - alpha/2, na.rm = TRUE),
      risk1_lo = quantile(risk1, alpha/2,     na.rm = TRUE),
      risk1_hi = quantile(risk1, 1 - alpha/2, na.rm = TRUE),
      rd_lo    = quantile(risk1 - .data[[risk0_var]], alpha/2,     na.rm = TRUE),
      rd_hi    = quantile(risk1 - .data[[risk0_var]], 1 - alpha/2, na.rm = TRUE),
      rr_lo    = quantile(risk1 / .data[[risk0_var]], alpha/2,     na.rm = TRUE),
      rr_hi    = quantile(risk1 / .data[[risk0_var]], 1 - alpha/2, na.rm = TRUE),
      .groups = "drop"
    )
  point.est |>
    mutate(risk0 = if (risk0_var == "risk0corr") risk0corr else risk0,
           rd = risk1 - risk0, rr = risk1 / risk0) |>
    dplyr::select(-dplyr::any_of(c("sim", "risk0corr"))) |>
    left_join(ci, by = "time_end")
}


# Additive results (point est + bootstrap) --------------------------------


## EQC additive: single bootstrap file
eqc.add.pt   <- readRDS(paste0(res_path, "eqc.add.pp.risk.pointest.rds"))
eqc.add.boot <- readRDS(paste0(res_path, "eqc.add.pp.boot.long.rds"))
eqc.add.ci   <- add.boot.ci(eqc.add.pt, eqc.add.boot)
saveRDS(eqc.add.ci, paste0(res_path, "eqc.add.pp.risks.ci.rds"))

## PCI additive: combine per-rep bootstrap files
pci_add_rep_files <- list.files(
  paste0(res_path, "pci_add_boot_reps"),
  pattern = "^pci_add_pp_boot_rep_\\d{3}\\.rds$",
  full.names = TRUE
)
pci.add.boot <- pci_add_rep_files |>
  lapply(readRDS) |>
  lapply(as.data.frame) |>
  bind_rows() |>
  as_tibble() |>
  mutate(across(everything(), as.numeric))
pci.add.pt <- readRDS(paste0(res_path, "pci.add.pp.risk.pointest.rds"))
pci.add.ci <- add.boot.ci(pci.add.pt, pci.add.boot)
saveRDS(pci.add.ci, paste0(res_path, "pci.add.pp.risks.ci.rds"))


# Multiplicative results (already have CIs) -------------------------------


eqc.mult.ci <- readRDS(paste0(res_path, "eqc.pp.risks.ci.rds"))
pci.mult.ci <- readRDS(paste0(res_path, "pci.pp.risks.ci.rds"))


# Overlay plot ------------------------------------------------------------

## solid = multiplicative (pooled logistic); dashed = additive (pooled linear)
## teal = no booster; orange = booster

overlay <- function(mult, add, title, ylim = c(0, 0.075)){
  mult <- mult |> arrange(time_end)
  add  <- add  |> arrange(time_end)
  col0 <- "#006663"; col1 <- "#FF6B1A"

  par(mar = c(5.1, 5.5, 4.1, 2.1))
  plot(NULL, xlim = range(c(0, mult$time_end)), ylim = ylim,
       xlab = "Weeks", ylab = "Risk", main = title,
       cex.axis = 1.4, cex.lab = 1.5, cex.main = 1.4, font.main = 1)
  grid()

  ## multiplicative ribbons + solid lines
  polygon(c(mult$time_end, rev(mult$time_end)), c(mult$risk0_lo, rev(mult$risk0_hi)),
          col = adjustcolor(col0, alpha.f = 0.18), border = NA)
  polygon(c(mult$time_end, rev(mult$time_end)), c(mult$risk1_lo, rev(mult$risk1_hi)),
          col = adjustcolor(col1, alpha.f = 0.18), border = NA)
  lines(c(0, mult$time_end), c(0, mult$risk0), col = col0, lty = 1, lwd = 2)
  lines(c(0, mult$time_end), c(0, mult$risk1), col = col1, lty = 1, lwd = 2)

  ## additive ribbons + dashed lines
  polygon(c(add$time_end, rev(add$time_end)), c(add$risk0_lo, rev(add$risk0_hi)),
          col = adjustcolor(col0, alpha.f = 0.10), border = NA)
  polygon(c(add$time_end, rev(add$time_end)), c(add$risk1_lo, rev(add$risk1_hi)),
          col = adjustcolor(col1, alpha.f = 0.10), border = NA)
  lines(c(0, add$time_end), c(0, add$risk0), col = col0, lty = 2, lwd = 2)
  lines(c(0, add$time_end), c(0, add$risk1), col = col1, lty = 2, lwd = 2)

  legend("topleft",
         legend = c("No Booster", "Booster", "Multiplicative (logistic)", "Additive (linear)"),
         col = c(col0, col1, "black", "black"),
         lty = c(1, 1, 1, 2), lwd = 2, cex = 1.1, bty = "n")
}

png(paste0(plot_path, "additive_vs_mult.cif.pp.png"), width = 2400, height = 4000, res = 300)
layout(matrix(1:2, nrow = 2))
overlay(eqc.mult.ci, eqc.add.ci, "Equi-confounding: additive vs multiplicative CIF")
overlay(pci.mult.ci, pci.add.ci, "Proximal inference: additive vs multiplicative CIF")
dev.off()
