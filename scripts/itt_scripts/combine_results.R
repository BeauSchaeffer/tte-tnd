##----- Kaiser Causal TTE-TND
##----- Combine Results and Draft Figures


# Packages ----------------------------------------------------------------


library(tidyverse)


# Data --------------------------------------------------------------------


res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results_itt.2/"

# STD Cox

  std.itt.cox.pointest <- readRDS(paste0(res_path,"std.itt.cox.pointest.rds"))
  std.itt.cox.nco.pointest <- readRDS(paste0(res_path,"std.cox.nco.itt.tidy.rds"))

# STD Pooled

  std.itt.risk.pointest <- readRDS(paste0(res_path, "std.itt.risk.pointest.rds"))
  std.itt.boot.long <- readRDS(paste0(res_path, "std.itt.boot.long.rds"))

# TND

  tnd.pointest <- readRDS(paste0(res_path,"tnd.itt.pointest.rds"))

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
  
  pci_rep_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results_weekmatch/pci_boot_reps/"
  
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
  
  # saveRDS(pci.itt.boot.long, file.path(res_path, "pci.itt.boot.long.rds")) # 2026-04-06
  
  
  
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


# STD ITT draft plot ------------------------------------------------------


std.itt.risks.ci <- pooled.boot.ci(point.est = std.itt.risk.pointest, boot.long = std.itt.boot.long)
# saveRDS(std.itt.risks.ci, paste0(res_path, "std.itt.risks.ci.rds")) # 2026-04-06
# png("figures_draft_wm/std.itt.risks.ci.plot.png", width = 2400, height=1800, res=300)
plot.risk.with.boot.ci(std.itt.risks.ci, title.main  = "Measured Covariate Adjustment Approach", title.sub = NULL)
# dev.off() # 2026-06-11


# EQC ITT draft plot ------------------------------------------------------


eqc.itt.HRs.ci <- eqc.cox.boot.ci(eqc.itt.cox.pointest, eqc.itt.cox.boot.long)
# saveRDS(eqc.itt.HRs.ci, paste0(res_path, "eqc.itt.HRs.ci.rds")) # 2026-04-06

eqc.itt.risks.ci <- pooled.boot.ci(point.est = eqc.itt.risk.pointest, boot.long = eqc.itt.boot.long)
# saveRDS(eqc.itt.risks.ci, paste0(res_path, "eqc.itt.risks.ci.rds")) # 2026-04-06

# png("figures_draft_wm/eqc.itt.risks.ci.plot.png", width = 2400, height=1800, res=300)
plot.risk.with.boot.ci(eqc.itt.risks.ci, title.main  = "Equi-confounding Approach", title.sub = NULL)
# dev.off() # 2026-06-11


# PCI ITT draft plot ------------------------------------------------------

pci.itt.HRs.ci <- pci.cox.boot.ci(pci.itt.cox.pointest, pci.itt.cox.boot.long)
# saveRDS(pci.itt.HRs.ci, paste0(res_path, "pci.itt.HRs.ci.rds")) # 2026-04-06

pci.itt.risks.ci <- pooled.boot.ci(point.est = pci.itt.risk.pointest, boot.long = pci.itt.boot.long)
# saveRDS(pci.itt.risks.ci, paste0(res_path, "pci.itt.risks.ci.rds")) # 2026-04-06

# png("figures_draft_wm/pci.itt.risks.ci.plot.png", width = 2400, height=1800, res=300)
plot.risk.with.boot.ci(pci.itt.risks.ci, title.main  = "Proximal inference Approach", title.sub = NULL)
# dev.off() # 2026-06-11


# STE and test behavior multipanel plot -----------------------------------

library(tidycmprsk)

data_Y3 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch/data_Y3_weekmatch.rds")

data_Y3 <- data_Y3 |> 
  mutate(Y3_itt_factor = case_when(
    Y3_itt_trunc==0 ~ "Censor",
    Y3_itt_trunc==1 ~ "Test Negative",
    Y3_itt_trunc==2 ~ "Test Positive"
  )) |> 
  mutate(Y3_itt_factor = factor(Y3_itt_factor, levels = c("Censor", "Test Negative", "Test Positive")),
         subclass=as.character(subclass))

eqc_Y3_cif_itt <- cuminc(
  Surv(Y3_itt_t_trunc, Y3_itt_factor) ~ treatment,
  data = data_Y3
)
# names(eqc_Y3_cif_itt$tidy)


# png("figures_draft_wm/multi.ste.testbehav.plot.png", width = 2400, height=4000, res=300)
par(mar = c(5.1, 4.1, 4.1, 2.1))
layout(matrix(1:2, nrow = 2))

## --- Panel 1: test pos risk curves under ste design  ---

plot.risk.with.boot.ci(
  std.itt.risks.ci,
  title.main = "Measured Covariate Adjustment Approach",
  title.sub  = NULL
)
mtext("A", side=3, adj=0, line=2, cex=1.5, font=1)

## --- Panel 2: test negative risk curves ---

par(mar = c(5.1, 5.5, 4.1, 2.1))
plot(NULL,
     xlim = range(c(0, eqc_Y3_cif_itt$tidy$time)),
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
tn0 <- eqc_Y3_cif_itt$tidy[eqc_Y3_cif_itt$tidy$outcome == "Test Negative" & eqc_Y3_cif_itt$tidy$strata == 0, ]
tn1 <- eqc_Y3_cif_itt$tidy[eqc_Y3_cif_itt$tidy$outcome == "Test Negative" & eqc_Y3_cif_itt$tidy$strata == 1, ]
# lines(c(eqc_Y3_cif_itt$tidy$time[eqc_Y3_cif_itt$tidy$outcome=="Test Positive" & eqc_Y3_cif_itt$tidy$strata==0]),
#       c(eqc_Y3_cif_itt$tidy$estimate[eqc_Y3_cif_itt$tidy$outcome=="Test Positive" & eqc_Y3_cif_itt$tidy$strata==0]),
#       col='#006663', lty=1, lwd=2)
# lines(c(eqc_Y3_cif_itt$tidy$time[eqc_Y3_cif_itt$tidy$outcome=="Test Positive" & eqc_Y3_cif_itt$tidy$strata==1]),
#       c(eqc_Y3_cif_itt$tidy$estimate[eqc_Y3_cif_itt$tidy$outcome=="Test Positive" & eqc_Y3_cif_itt$tidy$strata==1]),
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
lines(c(eqc_Y3_cif_itt$tidy$time[eqc_Y3_cif_itt$tidy$outcome=="Test Negative" & eqc_Y3_cif_itt$tidy$strata==0]),
      c(eqc_Y3_cif_itt$tidy$estimate[eqc_Y3_cif_itt$tidy$outcome=="Test Negative" & eqc_Y3_cif_itt$tidy$strata==0]),
      col='#006663', lty=2, lwd=2)
lines(c(eqc_Y3_cif_itt$tidy$time[eqc_Y3_cif_itt$tidy$outcome=="Test Negative" & eqc_Y3_cif_itt$tidy$strata==1]),
      c(eqc_Y3_cif_itt$tidy$estimate[eqc_Y3_cif_itt$tidy$outcome=="Test Negative" & eqc_Y3_cif_itt$tidy$strata==1]),
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
mtext("B", side=3, adj=0, line=2, cex=1.5, font=1)
# dev.off() # 2026-06-11






