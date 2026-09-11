##----- Kaiser Causal TTE-TND
##----- Simulation study: net-risk (Option A/B) vs competing-risk EQC de-biasing
##----- DRIVER: runs replications, summarizes bias / SD / RMSE, saves results
##-----
##----- Usage:  Rscript code/scripts/sim/sim_netrisk_run.R [N] [n_rep] [seed]
##-----
##----- Expected pattern (true b2A is a constant; frailty makes marginal hazards non-PH):
##-----   * de-biased b2A is recovered by (competing-risk OR net Option B) x (pooled OR cox_tv);
##-----   * net Option A is biased (survivor-conditioning mismatch across risk sets);
##-----   * scalar Cox is biased in every formulation (it collapses the time-varying
##-----     marginal log-HR, so the frailty attenuation does not cancel on subtraction);
##-----   * pooled CIF hits its formulation-specific oracle (net for A/B, crude for CR),
##-----     and risk0corr (de-biased) beats risk0 (naive) toward the truth.


here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) "code/scripts/sim")
source(file.path(here, "sim_netrisk_functions.R"))

args  <- commandArgs(trailingOnly = TRUE)
N     <- if (length(args) >= 1) as.integer(args[1]) else 5000L
n_rep <- if (length(args) >= 2) as.integer(args[2]) else 500L
seed  <- if (length(args) >= 3) as.integer(args[3]) else 20260818L

set.seed(seed)
message(sprintf("Simulation: N=%d, n_rep=%d, seed=%d", N, n_rep, seed))

b2A_all <- vector("list", n_rep)
cif_all <- vector("list", n_rep)

t0     <- Sys.time()
pb     <- txtProgressBar(min = 0, max = n_rep, style = 3)
n_fail <- 0L
for (r in seq_len(n_rep)) {
  out <- tryCatch(run_once(N), error = function(e) NULL)
  if (is.null(out)) {
    n_fail <- n_fail + 1L
  } else {
    b2A_all[[r]] <- cbind(rep = r, out$b2A)
    cif_all[[r]] <- cbind(rep = r, out$cif)
  }
  setTxtProgressBar(pb, r)
}
close(pb)
cat(sprintf("\nDone: %d reps (%d failed) in %.1f min\n", n_rep, n_fail,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

b2A_long <- data.table::rbindlist(b2A_all)
cif_long <- data.table::rbindlist(cif_all)

if (nrow(b2A_long) == 0) stop("All replications failed; check the DGP/estimator functions.")


# Summaries ---------------------------------------------------------------


## b2A: bias, empirical SD, RMSE by (formulation, estimator)
b2A_summary <- b2A_long[, .(
  n      = .N,
  truth  = mean(b2A_true),
  est    = mean(b2A),
  bias   = mean(b2A - b2A_true),
  sd     = sd(b2A),
  rmse   = sqrt(mean((b2A - b2A_true)^2))
), by = .(formulation, estimator)][order(formulation, estimator)]

## pooled CIF at landmark tau vs formulation-specific oracle
cif_summary <- cif_long[, .(
  n           = .N,
  risk1       = mean(risk1),      truth1 = mean(truth1),  bias1     = mean(risk1 - truth1),
  risk0_naive = mean(risk0),      truth0 = mean(truth0),  bias0_nv  = mean(risk0 - truth0),
  risk0_corr  = mean(risk0corr),                          bias0_cor = mean(risk0corr - truth0)
), by = .(formulation, estimator)][order(formulation, estimator)]

cat("\n================ b2A recovery (true = ", round(default_params$b2A, 3), ") ================\n", sep = "")
print(b2A_summary, digits = 3)

cat("\n================ CIF at week ", default_params$tau,
    " by estimator (naive vs de-biased risk0) ================\n", sep = "")
print(cif_summary, digits = 3)


# Save --------------------------------------------------------------------


out_dir <- here
saveRDS(list(b2A_long = b2A_long, cif_long = cif_long,
             b2A_summary = b2A_summary, cif_summary = cif_summary,
             params = default_params, N = N, n_rep = n_rep, seed = seed),
        file.path(out_dir, "sim_netrisk_results.rds"))
message("Saved: ", file.path(out_dir, "sim_netrisk_results.rds"))
