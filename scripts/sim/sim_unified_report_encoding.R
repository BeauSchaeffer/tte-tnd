#!/usr/bin/env Rscript
##----- Kaiser Causal TTE-TND
##----- UNIFIED Simulation Study: Report 3 - Risk-Set Encoding Robustness
##----- Demonstrates that results are (largely) invariant to the choice of risk-set
##----- encoding (competing-risks vs. net-risk Option A vs. net-risk Option B),
##----- using the pooled logistic estimator throughout.
##-----

library(data.table)

date <- "20260819"
results_file <- paste0("code/scripts/sim/sim_unified_b2A_", date, ".rds")

cat("Loading results from:", results_file, "\n")
results <- readRDS(results_file)

## Pooled logistic only, across all three encodings
fig_data <- results[estimator == "pooled"]
fig_data[, strategy := factor(strategy, levels = c("naive", "eqc", "pci"))]
fig_data[, encoding := factor(encoding, levels = c("CR", "netA", "netB"))]

fig_summary <- fig_data[, .(bias = mean(b2A - b2A_true)),
                         by = .(scenario, strategy, encoding)]

## LaTeX table: bias by scenario x strategy x encoding (pooled estimator)
tex_data <- fig_data[, .(
  bias = round(mean(b2A - b2A_true), 4),
  se   = round(sd(b2A), 4),
  rmse = round(sqrt(mean((b2A - b2A_true)^2)), 4)
), by = .(scenario, strategy, encoding)]

scen_ord <- c("eqc_holds", "eqc_violated")
enc_ord  <- c("CR", "netA", "netB")
enc_lab  <- c(CR = "Competing risk", netA = "Net-risk, Option A", netB = "Net-risk, Option B")
meth_lab <- c(naive = "Naive", eqc = "Equi-confounding", pci = "Proximal (PCI)")
meth_lab_short <- c(naive = "Naive", eqc = "EQC", pci = "PCI")
col_enc  <- c(CR = "#006663", netA = "#4C4CFF", netB = "#FF6B1A")

tex_data <- tex_data[order(match(scenario, scen_ord), match(strategy, c("naive", "eqc", "pci")),
                           match(encoding, enc_ord))]

cat("\n\nRisk-Set Encoding Robustness (pooled logistic, 200 reps each):\n")
print(tex_data)

tex_file <- "code/scripts/sim/tab_unified_encoding.tex"
sink(tex_file)
cat("\\begin{table}\n")
cat("\\centering\n")
cat("\\begin{tabular}{lllrrr}\n")
cat("\\toprule\n")
cat("Scenario & Strategy & Encoding & Bias & SE & RMSE \\\\\n")
cat("\\midrule\n")

prev_scenario <- ""; prev_strategy <- ""
for (i in seq_len(nrow(tex_data))) {
  row <- tex_data[i]
  scen_label <- ifelse(row$scenario == "eqc_holds", "Equi-confounding holds", "Equi-confounding violated")
  new_block  <- row$scenario != prev_scenario | row$strategy != prev_strategy
  if (new_block && i > 1) cat("\\addlinespace\n")
  cat(sprintf("%s & %s & %s & %.4f & %.4f & %.4f \\\\\n",
              ifelse(row$scenario != prev_scenario, scen_label, ""),
              ifelse(new_block, meth_lab[[as.character(row$strategy)]], ""),
              enc_lab[[as.character(row$encoding)]],
              row$bias, row$se, row$rmse))
  prev_scenario <- row$scenario; prev_strategy <- row$strategy
}

cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\caption{Robustness to the choice of risk-set encoding: bias, standard error, and RMSE of the causal log hazard-ratio estimate under three de-biasing strategies and three risk-set encodings (competing-risks, net-risk Option A, net-risk Option B), by equi-confounding assumption status. All estimates use the pooled discrete-time logistic model. Results based on 200 replications with $N = 4{,}000$ per replication.}\n")
cat("\\label{tab:unified_encoding}\n")
cat("\\end{table}\n")
sink()

cat("\nLaTeX table written to:", tex_file, "\n")

## Figure: Sampling distributions across encodings and strategies (base-R grouped boxplot)
strat_ord <- c("naive", "eqc", "pci")
pos <- as.vector(outer(0:2, (seq_along(strat_ord) - 1) * 4, "+")) + 1
grp_ord <- as.vector(outer(enc_ord, strat_ord, paste, sep = "."))

draw_fig <- function() {
  par(mfrow = c(1, 2), mar = c(4, 4.5, 2.5, 1))
  for (s in scen_ord) {
    d <- fig_data[scenario == s]
    d[, grp := factor(paste(encoding, strategy, sep = "."), levels = grp_ord)]
    boxplot(b2A ~ grp, data = d, at = pos, xaxt = "n", outline = FALSE,
            col = adjustcolor(col_enc[rep(enc_ord, length(strat_ord))], 0.6),
            xlab = "", ylab = expression(hat(beta)[2 * A]),
            main = c(eqc_holds = "Equi-confounding holds",
                     eqc_violated = "Equi-confounding violated")[[s]], cex.axis = 0.9,
            boxwex = 0.8)
    axis(1, at = pos[seq(2, length(pos), 3)], labels = meth_lab_short[strat_ord], cex.axis = 0.9)
    abline(h = -0.7, lty = 2, lwd = 2, col = "grey30")
    legend("topright", legend = enc_lab[enc_ord], fill = adjustcolor(col_enc[enc_ord], 0.6),
           bty = "n", cex = 0.8)
  }
}
png("code/scripts/sim/fig_unified_encoding.png", width = 2400, height = 1200, res = 240); draw_fig(); dev.off()
pdf("code/scripts/sim/fig_unified_encoding.pdf", width = 10, height = 5); draw_fig(); dev.off()

cat("Figures written to code/scripts/sim/fig_unified_encoding.{png,pdf}\n")

## Summary message
cat("\n\n===== RISK-SET ENCODING ROBUSTNESS SUMMARY =====\n")
for (sc in scen_ord) {
  cat("Scenario:", ifelse(sc == "eqc_holds", "Equi-confounding HOLDS", "Equi-confounding VIOLATED"), "\n")
  for (st in c("naive", "eqc", "pci")) {
    sub <- fig_summary[scenario == sc & strategy == st]
    cat(sprintf("  %s: CR=%.4f, netA=%.4f, netB=%.4f\n",
                st,
                sub[encoding == "CR"]$bias,
                sub[encoding == "netA"]$bias,
                sub[encoding == "netB"]$bias))
  }
}
cat("\n")
