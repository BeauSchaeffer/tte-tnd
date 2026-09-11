#!/usr/bin/env Rscript
##----- Kaiser Causal TTE-TND
##----- UNIFIED Simulation Study: CIF curve sampling distributions
##----- Generates CIF curves over time (weeks 1-53) for naive/eqc/pci strategies
##----- CR encoding, pooled logistic estimator, both scenarios
##----- Uses fewer replications (100) since only quantile bands are needed
##-----

library(data.table)

n_rep_curve <- 100
N <- 4000
seed0 <- 20260819

source("code/scripts/sim/sim_unified_functions.R")

all_curves <- list()
all_truth  <- list()

for (scenario in c("eqc_holds", "eqc_violated")) {
  cat("Scenario:", scenario, "\n")
  pb <- txtProgressBar(max = n_rep_curve, style = 3)

  for (i in seq_len(n_rep_curve)) {
    set.seed(seed0 + 10000 + i)
    g  <- gen_data(N, default_params, scenario); dd <- g$dd; negs <- g$negs
    p  <- g$p; tau <- p$tau

    P_cr <- pooled_cr(dd, tau)

    naive_p <- est_naive_pooled(P_cr$pos, tau)
    eqc_p   <- est_eqc_pooled(P_cr$pos, P_cr$neg, tau)
    pci_p   <- est_pci_pooled(P_cr$pos, P_cr$neg, tau)

    c_naive <- pooled_cif_curve(NULL, naive_p$m2, dd, tau, "naive")
    c_eqc   <- pooled_cif_curve(eqc_p$m1, eqc_p$m2, dd, tau, "eqc")
    c_pci   <- pooled_cif_curve(pci_p$m1, pci_p$m2, dd, tau, "pci")
    c_truth <- oracle_att_curve(dd, p, tau)

    all_curves[[length(all_curves) + 1]] <- rbind(
      cbind(c_naive, strategy = "naive", rep = i, scenario = scenario),
      cbind(c_eqc,   strategy = "eqc",   rep = i, scenario = scenario),
      cbind(c_pci,   strategy = "pci",   rep = i, scenario = scenario)
    )
    all_truth[[length(all_truth) + 1]] <- cbind(c_truth, rep = i, scenario = scenario)
    setTxtProgressBar(pb, i)
  }
  close(pb)
}

curves <- rbindlist(all_curves)
truth  <- rbindlist(all_truth)
saveRDS(list(curves = curves, truth = truth), "code/scripts/sim/sim_unified_cif_curves.rds")
cat("\nCurve data saved:", nrow(curves), "rows\n")

## Quantile bands (10th, 50th, 90th percentile) of risk0corr across replications
qbands <- curves[, .(
  q10 = quantile(risk0corr, 0.10),
  q50 = quantile(risk0corr, 0.50),
  q90 = quantile(risk0corr, 0.90)
), by = .(scenario, strategy, time)]

## Mean oracle truth curve across replications
truth_mean <- truth[, .(F0_true = mean(F0_true)), by = .(scenario, time)]


# Figure (base-R, matching color scheme of fig_unified_main.png) ----------


scen_ord <- c("eqc_holds", "eqc_violated")
scen_lab <- c(eqc_holds = "Equi-confounding holds", eqc_violated = "Equi-confounding violated")
meth_lab <- c(naive = "Naive", eqc = "EQC", pci = "PCI")
strat_ord <- c("naive", "eqc", "pci")
cols     <- c(naive = "#B22222", eqc = "#006663", pci = "#FF6B1A")

draw_fig <- function() {
  par(mfrow = c(1, 2), mar = c(4, 4.5, 2.5, 1))
  for (s in scen_ord) {
    qb <- qbands[scenario == s]
    tr <- truth_mean[scenario == s]
    ylim <- range(0, qb$q90, tr$F0_true)
    plot(NA, xlim = range(qb$time), ylim = ylim, xlab = "Week",
         ylab = "Cumulative incidence (untreated counterfactual, ATT)",
         main = scen_lab[[s]], cex.axis = 0.9)
    for (st in strat_ord) {
      d <- qb[strategy == st][order(time)]
      polygon(c(d$time, rev(d$time)), c(d$q10, rev(d$q90)),
              col = adjustcolor(cols[[st]], 0.2), border = NA)
    }
    for (st in strat_ord) {
      d <- qb[strategy == st][order(time)]
      lines(d$time, d$q50, col = cols[[st]], lwd = 2)
    }
    lines(tr$time[order(tr$time)], tr$F0_true[order(tr$time)], lty = 2, lwd = 2, col = "grey20")
    legend("topleft", legend = c(meth_lab[strat_ord], "True F0"),
           col = c(cols[strat_ord], "grey20"), lwd = 2, lty = c(1, 1, 1, 2),
           bty = "n", cex = 0.8)
  }
}
png("code/scripts/sim/fig_unified_cif_curves.png", width = 2400, height = 1200, res = 240); draw_fig(); dev.off()
pdf("code/scripts/sim/fig_unified_cif_curves.pdf", width = 10, height = 5); draw_fig(); dev.off()

cat("Figure written to code/scripts/sim/fig_unified_cif_curves.{png,pdf}\n")

