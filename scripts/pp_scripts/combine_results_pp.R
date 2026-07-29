##----- Kaiser Causal TTE-TND
##----- Combine Results and Draft Figures
##----- Per-protocol, no censoring weights


# Packages ----------------------------------------------------------------


library(tidyverse)


# Data --------------------------------------------------------------------


res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results_pp.3/"
plot_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/plots_pp.3/"

# STD Cox

  std.pp.cox.pointest <- readRDS(paste0(res_path,"std.pp.cox.pointest.rds"))
  std.pp.cox.nco.pointest <- readRDS(paste0(res_path,"std.cox.nco.pp.tidy.rds"))

# STD Pooled

  std.pp.risk.pointest <- readRDS(paste0(res_path, "std.pp.risk.pointest.rds"))
  std.pp.boot.long <- readRDS(paste0(res_path, "std.pp.boot.long.rds"))

# TND

  tnd.pp.pointest <- readRDS(paste0(res_path,"tnd.pp.pointest.rds"))

# EQC Cox

  eqc.pp.cox.pointest <- readRDS(paste0(res_path, "eqc.pp.cox.pointest.rds"))
  eqc.pp.cox.boot.long <- readRDS(paste0(res_path, "eqc.pp.cox.boot.long.rds"))

# EQC Pooled

  eqc.pp.risk.pointest <- readRDS(paste0(res_path, "eqc.pp.risk.pointest.rds"))
  eqc.pp.boot.long <- readRDS(paste0(res_path, "eqc.pp.boot.long.rds"))

# PCI Cox

  pci.pp.cox.pointest <- readRDS(paste0(res_path, "pci.pp.cox.pointest.rds"))
  pci.pp.cox.boot.long <- readRDS(paste0(res_path, "pci.pp.cox.boot.long.rds"))

# PCI Pooled
  pci.pp.risk.pointest <- readRDS(paste0(res_path, "pci.pp.risk.pointest.rds"))

  pci_rep_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results_pp.3/pci_boot_reps"
  
  pci_rep_files <- list.files(
    pci_rep_path,
    pattern = "^pci_pp_boot_rep_\\d{3}\\.rds$",
    full.names = TRUE
  )
  
  pci.pp.boot.long <- pci_rep_files |>
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
  
  saveRDS(pci.pp.boot.long, file.path(res_path, "pci.pp.boot.long.rds"))

  pci.pp.boot.long <- readRDS(paste0(res_path, "pci.pp.boot.long.rds"))


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
      
      # derived measures computed per bootstrap draw, then quantiled
      rd_lo   = quantile(risk1 - .data[[risk0_var]], probs = alpha/2, na.rm = TRUE),
      rd_hi   = quantile(risk1 - .data[[risk0_var]], probs = 1 - alpha/2, na.rm = TRUE),
      rr_lo   = quantile(risk1 / .data[[risk0_var]], probs = alpha/2, na.rm = TRUE),
      rr_hi   = quantile(risk1 / .data[[risk0_var]], probs = 1 - alpha/2, na.rm = TRUE),
      
      .groups = "drop"
    )
  
  # ensure point estimate uses the same risk0 definition
  point.est.use <- point.est |>
    mutate(risk0 = if (risk0_var == "risk0corr") risk0corr else risk0) |> 
    mutate(
      rd = risk1 - risk0,
      rr = risk1 / risk0
    )
  
  boot.point.ci <- boot.ci |>
    left_join(point.est.use, by = "time_end") |>
    select(-sim, -dplyr::any_of("risk0corr")) |>
    relocate(
      time_end,
      risk0, risk0_lo, risk0_hi,
      risk1, risk1_lo, risk1_hi,
      rd, rd_lo, rd_hi,
      rr, rr_lo, rr_hi
    )
  
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
    treatHR = unname(eqc.pp.cox.pointest["treatHR"]),
    fluvaxHR = unname(eqc.pp.cox.pointest["fluvaxHR"])
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
       cex.main = cex.main,
       font.main = 1)
  
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


# STD PP draft plot ------------------------------------------------------


std.pp.risks.ci <- pooled.boot.ci(point.est = std.pp.risk.pointest, boot.long = std.pp.boot.long)
saveRDS(std.pp.risks.ci, paste0(res_path, "std.pp.risks.ci.rds"))

png(paste0(plot_path,"std.pp.risks.ci.plot.png"), width = 2400, height=1800, res=300)
plot.risk.with.boot.ci(std.pp.risks.ci, title.main = "Measured Covariate Adjustment Approach", title.sub = NULL)
dev.off()


# EQC PP draft plot ------------------------------------------------------


eqc.pp.HRs.ci <- eqc.cox.boot.ci(eqc.pp.cox.pointest, eqc.pp.cox.boot.long)
saveRDS(eqc.pp.HRs.ci, paste0(res_path, "eqc.pp.HRs.ci.rds"))

eqc.pp.risks.ci <- pooled.boot.ci(point.est = eqc.pp.risk.pointest, boot.long = eqc.pp.boot.long)
saveRDS(eqc.pp.risks.ci, paste0(res_path, "eqc.pp.risks.ci.rds"))

png(paste0(plot_path,"eqc.pp.risks.ci.plot.png"), width = 2400, height=1800, res=300)
plot.risk.with.boot.ci(eqc.pp.risks.ci, title.main  = "Equi-confounding Approach", title.sub = NULL)
dev.off()


# PCI PP draft plot ------------------------------------------------------

pci.pp.HRs.ci <- pci.cox.boot.ci(pci.pp.cox.pointest, pci.pp.cox.boot.long)
saveRDS(pci.pp.HRs.ci, paste0(res_path, "pci.pp.HRs.ci.rds"))

pci.pp.risks.ci <- pooled.boot.ci(point.est = pci.pp.risk.pointest, boot.long = pci.pp.boot.long)
saveRDS(pci.pp.risks.ci, paste0(res_path, "pci.pp.risks.ci.rds"))

png(paste0(plot_path,"pci.pp.risks.ci.plot.png"), width = 2400, height=1800, res=300)
plot.risk.with.boot.ci(pci.pp.risks.ci, title.main  = "Proximal inference Approach", title.sub = NULL)
dev.off()


# STE and test behavior multipanel plot -----------------------------------


# Test behavior draft plot ------------------------------------------------


library(tidycmprsk)

data_Y3 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch.3/data_Y3_weekmatch.rds")

data_Y3 <- data_Y3 |> 
  mutate(Y3_pp_factor = case_when(
    Y3_pp_trunc==0 ~ "Censor",
    Y3_pp_trunc==1 ~ "Test Negative",
    Y3_pp_trunc==2 ~ "Test Positive"
  )) |> 
  mutate(Y3_pp_factor = factor(Y3_pp_factor, levels = c("Censor", "Test Negative", "Test Positive")),
         subclass=as.character(subclass))

eqc_Y3_cif_pp <- cuminc(
  Surv(Y3_pp_t_trunc, Y3_pp_factor) ~ treatment,
  data = data_Y3
)

# png("figures_draft_wm/multi.ste.testbehav.plot.png", width = 2400, height=4000, res=300)
png(paste0(plot_path,"testbehav.plot.png"), width = 2400, height=1800, res=300)
# par(mar = c(5.1, 4.1, 4.1, 2.1))
# layout(matrix(1:2, nrow = 2))

# ## --- Panel 1: test pos risk curves under ste design  ---
# 
# plot.risk.with.boot.ci(
#   std.pp.risks.ci,
#   title.main = "Measured Covariate Adjustment Approach",
#   title.sub  = NULL
# )
# mtext("A", side=3, adj=0, line=2, cex=1.5, font=1)

## --- Panel 2: test negative risk curves ---

par(mar = c(5.1, 5.5, 4.1, 2.1))
plot(NULL,
     xlim = range(c(0, eqc_Y3_cif_pp$tidy$time)),
     ylim = range(c(0, 0.15)),
     xlab="Weeks",
     ylab="Risk",
     main="Health-seeking Behavior",
     cex.axis = 1.5,
     cex.lab = 1.5,
     cex.main=1.4,
     font.main = 1
)
# mtext("Testing Behavior", side = 3, line = 0.5, font = 3, cex=1.2)
grid()
tn0 <- eqc_Y3_cif_pp$tidy[eqc_Y3_cif_pp$tidy$outcome == "Test Negative" & eqc_Y3_cif_pp$tidy$strata == 0, ]
tn1 <- eqc_Y3_cif_pp$tidy[eqc_Y3_cif_pp$tidy$outcome == "Test Negative" & eqc_Y3_cif_pp$tidy$strata == 1, ]
# lines(c(eqc_Y3_cif_pp$tidy$time[eqc_Y3_cif_pp$tidy$outcome=="Test Positive" & eqc_Y3_cif_pp$tidy$strata==0]),
#       c(eqc_Y3_cif_pp$tidy$estimate[eqc_Y3_cif_pp$tidy$outcome=="Test Positive" & eqc_Y3_cif_pp$tidy$strata==0]),
#       col='#006663', lty=1, lwd=2)
# lines(c(eqc_Y3_cif_pp$tidy$time[eqc_Y3_cif_pp$tidy$outcome=="Test Positive" & eqc_Y3_cif_pp$tidy$strata==1]),
#       c(eqc_Y3_cif_pp$tidy$estimate[eqc_Y3_cif_pp$tidy$outcome=="Test Positive" & eqc_Y3_cif_pp$tidy$strata==1]),
#       col='#FF6B1A', lty=1, lwd=2)
polygon(
  x = c(tn0$time, rev(tn0$time)),
  y = c(tn0$conf.low, rev(tn0$conf.high)),
  col = adjustcolor('#006663', alpha.f = 0.25),
  border = NA
)
polygon(
  x = c(tn1$time, rev(tn1$time)),
  y = c(tn1$conf.low, rev(tn1$conf.high)),
  col = adjustcolor('#FF6B1A', alpha.f = 0.25),
  border = NA
)
lines(c(eqc_Y3_cif_pp$tidy$time[eqc_Y3_cif_pp$tidy$outcome=="Test Negative" & eqc_Y3_cif_pp$tidy$strata==0]),
      c(eqc_Y3_cif_pp$tidy$estimate[eqc_Y3_cif_pp$tidy$outcome=="Test Negative" & eqc_Y3_cif_pp$tidy$strata==0]),
      col='#006663', lty=2, lwd=2)
lines(c(eqc_Y3_cif_pp$tidy$time[eqc_Y3_cif_pp$tidy$outcome=="Test Negative" & eqc_Y3_cif_pp$tidy$strata==1]),
      c(eqc_Y3_cif_pp$tidy$estimate[eqc_Y3_cif_pp$tidy$outcome=="Test Negative" & eqc_Y3_cif_pp$tidy$strata==1]),
      col='#FF6B1A', lty=2, lwd=2)
# legend("topleft",
#        legend = c("No Booster", "Booster"),
#        col = c('#006663', '#FF6B1A'),
#        lty = 1, lwd = 4, cex=1.2,
#        bty = "n")
# legend("topright",
#        legend = c("Test Positive", "Test Negative"),
#        col = 'black',
#        lty = c(1,2), lwd = 2, cex=1.2,
#        bty = "n")
# mtext("B", side=3, adj=0, line=2, cex=1.5, font=1)
dev.off()


# Forest plot -------------------------------------------------------------


forest.data <- bind_rows(
  # 1. Treatment HR on primary outcome
  std.pp.cox.pointest |>
    filter(term == "treatment") |>
    transmute(label = "HR Booster\nTest Positive",
              hr = estimate, lo = conf.low, hi = conf.high,
              group = "effect"),
  
  # 2. NCE (flu_vax) HR on primary outcome
  std.pp.cox.pointest |>
    filter(term == "flu_vax") |>
    transmute(label = "HR Flu Vax\nTest Negative",
              hr = estimate, lo = conf.low, hi = conf.high,
              group = "nc"),
  
  # 3. Treatment HR on NCO
  std.pp.cox.nco.pointest |>
    filter(term == "treatment") |>
    transmute(label = "HR Booster\nTest Negative",
              hr = estimate, lo = conf.low, hi = conf.high,
              group = "nc")
) |>
  mutate(y = rev(seq_along(label)),
         # col = if_else(group == "effect", "#FF6B1A", "#006663")
         col = "black")

# png(paste0(plot_path, "forest.plot.png"), width = 3200, height = 1800, res = 300)
# 
# par(mar = c(5.1, 13, 5.1, 2.1))
# 
# xlim <- range(c(1, forest.data$lo, forest.data$hi))
# xlim <- xlim + c(-0.15, 0.15) * diff(xlim)   # wider padding, no more edge clipping
# 
# plot(NULL,
#      xlim = xlim,
#      ylim = c(0.5, nrow(forest.data) + 0.7),  # extra headroom for text above top row
#      yaxt = "n",
#      xlab = "Hazard Ratio",
#      ylab = "",
#      main = "Main Effect and Negative Controls Estimates",
#      cex.axis = 1.4,
#      cex.lab  = 1.5,
#      cex.main = 1.4,
#      font.main = 1)
# 
# axis(2, at = forest.data$y, labels = forest.data$label, las = 1, cex.axis = 1.3)
# 
# abline(v = 1, lty = 2, col = "grey50")
# 
# segments(x0 = forest.data$lo, x1 = forest.data$hi,
#          y0 = forest.data$y, y1 = forest.data$y,
#          lwd = 2, col = forest.data$col)
# 
# points(forest.data$hr, forest.data$y, pch = 15, cex = 1.8, col = forest.data$col)
# 
# text(forest.data$hr, forest.data$y + 0.25,
#      labels = sprintf("%.2f (%.2f, %.2f)", forest.data$hr, forest.data$lo, forest.data$hi),
#      cex = 1.0, xpd = NA)
# 
# dev.off()


# STD and Forest Multipanel -----------------------------------------------

panel_label <- function(label, x_ndc = 0.02, y_frac = 0.97, ...) {
  usr <- par("usr")
  text(x = grconvertX(x_ndc, from = "ndc", to = "user"),
       y = usr[3] + y_frac * diff(usr[3:4]),
       labels = label, xpd = NA, adj = c(0, 0), ...)
}

png(paste0(plot_path, "std.pp.forest.multi.png"), width = 2400, height = 4000, res = 300)

layout(matrix(1:2, nrow = 2))

# --- Panel A: STD PP risk curves ---

par(mar = c(5.1, 5.5, 4.1, 2.1))
plot.risk.with.boot.ci(
  std.pp.risks.ci,
  title.main = "Measured Covariate Adjustment Approach",
  title.sub  = NULL
)
panel_label("A", cex = 1.5, font = 2)

## --- Panel B: Forest plot ---

par(mar = c(5.1, 13, 5.1, 2.1))
xlim <- range(c(1, forest.data$lo, forest.data$hi))
xlim <- xlim + c(-0.15, 0.15) * diff(xlim)

plot(NULL,
     xlim = xlim,
     ylim = c(0.5, nrow(forest.data) + 0.7),
     yaxt = "n",
     xlab = "Hazard Ratio",
     ylab = "",
     main = "Main Effect and Negative Controls Estimates",
     cex.axis = 1.4,
     cex.lab  = 1.5,
     cex.main = 1.4,
     font.main = 1)

axis(2, at = forest.data$y, labels = forest.data$label, las = 1, cex.axis = 1.3)

abline(v = 1, lty = 2, col = "grey50")

segments(x0 = forest.data$lo, x1 = forest.data$hi,
         y0 = forest.data$y, y1 = forest.data$y,
         lwd = 2, col = forest.data$col)

points(forest.data$hr, forest.data$y, pch = 15, cex = 1.8, col = forest.data$col)

text(forest.data$hr, forest.data$y + 0.25,
     labels = sprintf("%.2f (%.2f, %.2f)", forest.data$hr, forest.data$lo, forest.data$hi),
     cex = 1.0, xpd = NA)

panel_label("B", cex = 1.5, font = 2)

dev.off()





