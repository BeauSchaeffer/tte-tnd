#!/usr/bin/env Rscript
##----- Kaiser Causal TTE-TND
##----- UNIFIED Simulation Study: EQC vs PCI under assumption violation
##----- Runner script for full 2 scenarios × 2 encodings × 3 strategies × 3 estimators
##-----

library(data.table)

## Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
N      <- as.integer(args[1] %||% 10000)  # sample size per replication
n_rep  <- as.integer(args[2] %||% 500)    # number of replications
date   <- args[3] %||% format(Sys.Date(), "%Y%m%d")  # for results file naming
seed0  <- 20260819

## Source functions
source("code/scripts/sim/sim_unified_functions.R")

cat("\n==== Unified EQC vs PCI Simulation Study ====\n")
cat("Sample size per replication: N =", N, "\n")
cat("Number of replications: n_rep =", n_rep, "\n")
cat("Date tag:", date, "\n\n")

## Run simulations
all_b2A <- list()
all_cif <- list()
pb <- txtProgressBar(max = n_rep * 2, style = 3)  # 2 scenarios

for (scenario in c("eqc_holds", "eqc_violated")) {
  cat("\nScenario:", scenario, "\n")
  
  for (i in seq_len(n_rep)) {
    result <- run_once(N = N, scenario = scenario, seed = seed0 + i, cox_encodings = "CR")
    all_b2A[[length(all_b2A) + 1]] <- result$b2A
    all_cif[[length(all_cif) + 1]] <- result$cif
    
    if (i %% max(1, n_rep %/% 10) == 0) {
      setTxtProgressBar(pb, which(c("eqc_holds", "eqc_violated") == scenario) * n_rep - n_rep + i)
    }
  }
}
close(pb)

## Combine results
results_b2A <- rbindlist(all_b2A)
results_cif <- rbindlist(all_cif)

cat("\n\nSample sizes:\n")
cat("  b2A results:", nrow(results_b2A), "rows\n")
cat("  CIF results:", nrow(results_cif), "rows\n")

## Summary statistics by scenario x encoding x strategy x estimator
summary_b2A <- results_b2A[, .(
  bias  = mean(b2A - b2A_true),
  se    = sd(b2A),
  rmse  = sqrt(mean((b2A - b2A_true)^2))
), by = .(scenario, encoding, strategy, estimator)]

cat("\n\nSummary statistics (log-HR scale, bias/se/rmse):\n")
print(summary_b2A)

## Save results
saveRDS(results_b2A, paste0("code/scripts/sim/sim_unified_b2A_", date, ".rds"))
saveRDS(results_cif, paste0("code/scripts/sim/sim_unified_cif_", date, ".rds"))
saveRDS(summary_b2A, paste0("code/scripts/sim/sim_unified_summary_", date, ".rds"))

cat("\n\nResults saved:\n")
cat("  code/scripts/sim/sim_unified_b2A_", date, ".rds\n", sep="")
cat("  code/scripts/sim/sim_unified_cif_", date, ".rds\n", sep="")
cat("  code/scripts/sim/sim_unified_summary_", date, ".rds\n", sep="")

cat("\n\nDone!\n")
