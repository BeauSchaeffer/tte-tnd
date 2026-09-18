##----- Kaiser Causal TTE-TND
##----- Combine Results and Draft Figures -- NET RISK, Option B
##----- Per-protocol, no censoring weights
##----- Last updated 2026-09-16

# Packages ----------------------------------------------------------------


library(tidyverse)


# Data --------------------------------------------------------------------


res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results_pp.5/"
plot_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/plots_pp.5/"

dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)

# EQC Cox

  netB.eqc.pp.cox.pointest <- readRDS(paste0(res_path, "netB.eqc.pp.cox.pointest.rds"))
  netB.eqc.pp.cox.boot.long <- readRDS(paste0(res_path, "netB.eqc.pp.cox.boot.long.rds"))

# EQC Pooled

  netB.eqc.pp.risk.pointest <- readRDS(paste0(res_path, "netB.eqc.pp.risk.pointest.rds"))

  netB_eqc_rep_path <- paste0(res_path, "netB_eqc_boot_reps")

  netB_eqc_rep_files <- list.files(
    netB_eqc_rep_path,
    pattern = "^netB_eqc_pp_boot_rep_\\d{3}\\.rds$",
    full.names = TRUE
  )

  message("EQC bootstrap replicates found: ", length(netB_eqc_rep_files))

  netB.eqc.pp.boot.long <- netB_eqc_rep_files |>
    lapply(readRDS) |>
    lapply(\(m) as.data.frame(m)) |>
    bind_rows() |>
    as_tibble() |>
    mutate(
      sim       = as.integer(sim),
      time_end  = as.integer(time_end),
      risk0     = as.numeric(risk0),
      risk0corr = as.numeric(risk0corr),
      risk1     = as.numeric(risk1)
    ) |>
    arrange(sim, time_end)

  saveRDS(netB.eqc.pp.boot.long, file.path(res_path, "netB.eqc.pp.boot.long.rds"))

  netB.eqc.pp.boot.long <- readRDS(paste0(res_path, "netB.eqc.pp.boot.long.rds"))

# PCI Cox

  netB.pci.pp.cox.pointest <- readRDS(paste0(res_path, "netB.pci.pp.cox.pointest.rds"))
  netB.pci.pp.cox.boot.long <- readRDS(paste0(res_path, "netB.pci.pp.cox.boot.long.rds"))

# PCI Pooled

  netB.pci.pp.risk.pointest <- readRDS(paste0(res_path, "netB.pci.pp.risk.pointest.rds"))

  netB_pci_rep_path <- paste0(res_path, "netB_pci_boot_reps")

  netB_pci_rep_files <- list.files(
    netB_pci_rep_path,
    pattern = "^netB_pci_pp_boot_rep_\\d{3}\\.rds$",
    full.names = TRUE
  )

  message("PCI bootstrap replicates found: ", length(netB_pci_rep_files))

  netB.pci.pp.boot.long <- netB_pci_rep_files |>
    lapply(readRDS) |>
    lapply(\(m) as.data.frame(m)) |>
    bind_rows() |>
    as_tibble() |>
    mutate(
      sim      = as.integer(sim),
      time_end = as.integer(time_end),
      risk0    = as.numeric(risk0),
      risk1    = as.numeric(risk1)
    ) |>
    arrange(sim, time_end)

  saveRDS(netB.pci.pp.boot.long, file.path(res_path, "netB.pci.pp.boot.long.rds")) # run again after boots complete

  netB.pci.pp.boot.long <- readRDS(paste0(res_path, "netB.pci.pp.boot.long.rds"))


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
    treatHR = unname(point.est["treatHR"]),
    fluvaxHR = unname(point.est["fluvaxHR"])
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
                                   lty.lines = 3,
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
  
  # point estimate lines (dashed = net risk)
  lines(c(0, risks.and.cis$time_end),
        c(0, risks.and.cis$risk0),
        col = col0, lty = lty.lines, lwd = lwd.lines)
  
  lines(c(0, risks.and.cis$time_end),
        c(0, risks.and.cis$risk1),
        col = col1, lty = lty.lines, lwd = lwd.lines)
  
  legend(legend_pos,
         legend = c("No Booster", "Booster"),
         col = c(col0, col1),
         lty = lty.lines, lwd = lwd.legend, cex = legend.cex,
         bty = "n")
  
  invisible(risks.and.cis)
}


# EQC net risk draft plot ------------------------------------------------


netB.eqc.pp.HRs.ci <- eqc.cox.boot.ci(netB.eqc.pp.cox.pointest, netB.eqc.pp.cox.boot.long)
saveRDS(netB.eqc.pp.HRs.ci, paste0(res_path, "netB.eqc.pp.HRs.ci.rds"))

netB.eqc.pp.risks.ci <- pooled.boot.ci(point.est = netB.eqc.pp.risk.pointest,
                                       boot.long = netB.eqc.pp.boot.long)
saveRDS(netB.eqc.pp.risks.ci, paste0(res_path, "netB.eqc.pp.risks.ci.rds"))

png(paste0(plot_path,"netB.eqc.pp.risks.ci.plot.png"), width = 2400, height=1800, res=300)
plot.risk.with.boot.ci(netB.eqc.pp.risks.ci,
                       title.main = "Equi-confounding Approach",
                       title.sub  = "Net risk")
dev.off()


# PCI net risk draft plot ------------------------------------------------


netB.pci.pp.HRs.ci <- pci.cox.boot.ci(netB.pci.pp.cox.pointest, netB.pci.pp.cox.boot.long)
saveRDS(netB.pci.pp.HRs.ci, paste0(res_path, "netB.pci.pp.HRs.ci.rds"))

netB.pci.pp.risks.ci <- pooled.boot.ci(point.est = netB.pci.pp.risk.pointest,
                                       boot.long = netB.pci.pp.boot.long)
saveRDS(netB.pci.pp.risks.ci, paste0(res_path, "netB.pci.pp.risks.ci.rds"))

png(paste0(plot_path,"netB.pci.pp.risks.ci.plot.png"), width = 2400, height=1800, res=300)
plot.risk.with.boot.ci(netB.pci.pp.risks.ci,
                       title.main = "Proximal Inference Approach",
                       title.sub  = "Net risk")
dev.off()


# EQC and PCI net risk multipanel -----------------------------------------


panel_label <- function(label, x_ndc = 0.02, y_frac = 0.97, ...) {
  usr <- par("usr")
  text(x = grconvertX(x_ndc, from = "ndc", to = "user"),
       y = usr[3] + y_frac * diff(usr[3:4]),
       labels = label, xpd = NA, adj = c(0, 0), ...)
}

png(paste0(plot_path, "netB.pp.multi.png"), width = 2400, height = 4000, res = 300)

layout(matrix(1:2, nrow = 2))

## --- Panel A: EQC ---

plot.risk.with.boot.ci(
  netB.eqc.pp.risks.ci,
  title.main = "Equi-confounding Approach",
  title.sub  = "Net risk"
)
panel_label("A", cex = 1.5, font = 2)

## --- Panel B: PCI ---

plot.risk.with.boot.ci(
  netB.pci.pp.risks.ci,
  title.main = "Proximal Inference Approach",
  title.sub  = "Net risk"
)
panel_label("B", cex = 1.5, font = 2)

dev.off()


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

# Net risk vs competing risk overlay plotting function ---------------------

### Both sets of curves come from the data_weekmatch.3 cohort: identical
### matched pairs, identical index weeks. The two differ only in the
### definition of the risk set -- under competing events a negative test
### removes an individual from follow-up, under net risk only a positive
### does. This is the comparison the simulation study characterises.

plot.risk.overlay.boot.ci <- function(cr.risks.and.cis,
                                      net.risks.and.cis,
                                      title.main = "Risk Curves",
                                      title.sub  = "Net risk vs competing risk",
                                      xlab = "Weeks",
                                      ylab = "Risk",
                                      ylim = c(0, 0.075),
                                      col0 = "#006663",
                                      col1 = "#FF6B1A",
                                      cr.ribbon.alpha  = 0.18,
                                      net.ribbon.alpha = 0.10,
                                      lwd.lines = 2,
                                      lwd.legend = 4,
                                      legend_pos = "topleft",
                                      add.grid = TRUE,
                                      mar = c(5.1, 5.5, 4.1, 2.1),
                                      cex.axis = 1.5,
                                      cex.lab  = 1.5,
                                      cex.main = 1.4,
                                      cex.sub  = 1.2,
                                      legend.cex = 1.1) {
  
  req_cols <- c("time_end",
                "risk0", "risk0_lo", "risk0_hi",
                "risk1", "risk1_lo", "risk1_hi")
  
  missing.cr <- setdiff(req_cols, names(cr.risks.and.cis))
  if (length(missing.cr) > 0) {
    stop("`cr.risks.and.cis` is missing columns: ", paste(missing.cr, collapse = ", "))
  }
  
  missing.net <- setdiff(req_cols, names(net.risks.and.cis))
  if (length(missing.net) > 0) {
    stop("`net.risks.and.cis` is missing columns: ", paste(missing.net, collapse = ", "))
  }
  
  cr.risks.and.cis  <- cr.risks.and.cis  |> arrange(time_end)
  net.risks.and.cis <- net.risks.and.cis |> arrange(time_end)
  
  xlim <- range(c(0, cr.risks.and.cis$time_end, net.risks.and.cis$time_end), na.rm = TRUE)
  
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
  
  # competing risk CI ribbons
  polygon(
    x = c(cr.risks.and.cis$time_end, rev(cr.risks.and.cis$time_end)),
    y = c(cr.risks.and.cis$risk0_lo, rev(cr.risks.and.cis$risk0_hi)),
    col = adjustcolor(col0, alpha.f = cr.ribbon.alpha),
    border = NA
  )
  
  polygon(
    x = c(cr.risks.and.cis$time_end, rev(cr.risks.and.cis$time_end)),
    y = c(cr.risks.and.cis$risk1_lo, rev(cr.risks.and.cis$risk1_hi)),
    col = adjustcolor(col1, alpha.f = cr.ribbon.alpha),
    border = NA
  )
  
  # net risk CI ribbons
  polygon(
    x = c(net.risks.and.cis$time_end, rev(net.risks.and.cis$time_end)),
    y = c(net.risks.and.cis$risk0_lo, rev(net.risks.and.cis$risk0_hi)),
    col = adjustcolor(col0, alpha.f = net.ribbon.alpha),
    border = NA
  )
  
  polygon(
    x = c(net.risks.and.cis$time_end, rev(net.risks.and.cis$time_end)),
    y = c(net.risks.and.cis$risk1_lo, rev(net.risks.and.cis$risk1_hi)),
    col = adjustcolor(col1, alpha.f = net.ribbon.alpha),
    border = NA
  )
  
  # competing risk point estimate lines (solid)
  lines(c(0, cr.risks.and.cis$time_end),
        c(0, cr.risks.and.cis$risk0),
        col = col0, lty = 1, lwd = lwd.lines)
  
  lines(c(0, cr.risks.and.cis$time_end),
        c(0, cr.risks.and.cis$risk1),
        col = col1, lty = 1, lwd = lwd.lines)
  
  # net risk point estimate lines (dashed)
  lines(c(0, net.risks.and.cis$time_end),
        c(0, net.risks.and.cis$risk0),
        col = col0, lty = 3, lwd = lwd.lines)
  
  lines(c(0, net.risks.and.cis$time_end),
        c(0, net.risks.and.cis$risk1),
        col = col1, lty = 3, lwd = lwd.lines)
  
  legend(legend_pos,
         legend = c("No Booster", "Booster", "Competing risk", "Net risk"),
         col = c(col0, col1, "black", "black"),
         lty = c(1, 1, 1, 3),
         lwd = c(lwd.legend, lwd.legend, lwd.lines, lwd.lines),
         cex = legend.cex,
         bty = "n")
  
  invisible(list(cr = cr.risks.and.cis, net = net.risks.and.cis))
}


# Net risk vs competing risk overlays --------------------------------------


### CR curves from the .3 cohort, produced by combine_results_pp.R
eqc.pp.risks.ci <- readRDS(paste0(res_path, "eqc.pp.risks.ci.rds"))
pci.pp.risks.ci <- readRDS(paste0(res_path, "pci.pp.risks.ci.rds"))

png(paste0(plot_path, "netB.eqc.vs.cr.pp.risks.ci.plot.png"), width = 2400, height = 1800, res = 300)
plot.risk.overlay.boot.ci(
  cr.risks.and.cis  = eqc.pp.risks.ci,
  net.risks.and.cis = netB.eqc.pp.risks.ci,
  title.main = "Equi-confounding Approach",
  title.sub  = "Net risk vs competing risk"
)
dev.off()

png(paste0(plot_path, "netB.pci.vs.cr.pp.risks.ci.plot.png"), width = 2400, height = 1800, res = 300)
plot.risk.overlay.boot.ci(
  cr.risks.and.cis  = pci.pp.risks.ci,
  net.risks.and.cis = netB.pci.pp.risks.ci,
  title.main = "Proximal Inference Approach",
  title.sub  = "Net risk vs competing risk"
)
dev.off()


# Net risk vs competing risk multipanel ------------------------------------


png(paste0(plot_path, "netB.vs.cr.pp.multi.png"), width = 2400, height = 4000, res = 300)

layout(matrix(1:2, nrow = 2))

## --- Panel A: EQC ---

plot.risk.overlay.boot.ci(
  cr.risks.and.cis  = eqc.pp.risks.ci,
  net.risks.and.cis = netB.eqc.pp.risks.ci,
  title.main = "Equi-confounding Approach",
  title.sub  = "Net risk vs competing risk"
)
panel_label("A", cex = 1.5, font = 2)

## --- Panel B: PCI ---

plot.risk.overlay.boot.ci(
  cr.risks.and.cis  = pci.pp.risks.ci,
  net.risks.and.cis = netB.pci.pp.risks.ci,
  title.main = "Proximal Inference Approach",
  title.sub  = "Net risk vs competing risk"
)
panel_label("B", cex = 1.5, font = 2)

dev.off()


# Test behavior draft plot ------------------------------------------------


### Recurrent negative testing on the net risk risk set {T2 > t}. Unlike the
### competing risk version this is NOT a cuminc CIF -- negatives do not remove
### anyone from the risk set, so the natural summary is the cumulative mean
### number of negative tests per person by arm.

data_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch.3/"

data_Y2  <- read_rds(paste0(data_path, "data_Y2_weekmatch.rds"))
neg_hist <- read_rds(paste0(data_path, "neg_hist_weekmatch.rds"))

netB_neg_cum <- neg_hist |>
  inner_join(data_Y2 |> dplyr::select(fake_mrn, treatment, cap = Y2_pp_t_trunc),
             by = "fake_mrn") |>
  filter(neg_t <= cap) |>
  distinct(fake_mrn, treatment, time_end = neg_t) |>
  count(treatment, time_end, name = "n_neg") |>
  complete(treatment, time_end = 1:52, fill = list(n_neg = 0)) |>
  arrange(treatment, time_end) |>
  group_by(treatment) |>
  mutate(cum_neg = cumsum(n_neg)) |>
  ungroup() |>
  left_join(data_Y2 |> count(treatment, name = "n_arm"), by = "treatment") |>
  mutate(cum_neg_per_person = cum_neg / n_arm)

png(paste0(plot_path,"netB.testbehav.plot.png"), width = 2400, height=1800, res=300)

par(mar = c(5.1, 5.5, 4.1, 2.1))
plot(NULL,
     xlim = range(c(0, netB_neg_cum$time_end)),
     ylim = range(c(0, netB_neg_cum$cum_neg_per_person)),
     xlab = "Weeks",
     ylab = "Cumulative negative tests per person",
     main = "Health-seeking Behavior",
     cex.axis = 1.5,
     cex.lab = 1.5,
     cex.main = 1.4,
     font.main = 1
)
mtext("Recurrent testing on the net risk risk set", side = 3, line = 0.5, font = 3, cex = 1.2)
grid()

tn0 <- netB_neg_cum[netB_neg_cum$treatment == 0, ]
tn1 <- netB_neg_cum[netB_neg_cum$treatment == 1, ]

lines(c(0, tn0$time_end), c(0, tn0$cum_neg_per_person),
      col = '#006663', lty = 3, lwd = 2)
lines(c(0, tn1$time_end), c(0, tn1$cum_neg_per_person),
      col = '#FF6B1A', lty = 3, lwd = 2)

legend("topleft",
       legend = c("No Booster", "Booster"),
       col = c('#006663', '#FF6B1A'),
       lty = 3, lwd = 4, cex = 1.2,
       bty = "n")

dev.off()

