##----- Kaiser Causal TTE-TND
##----- Unified Simulation Study: EQC vs PCI under assumption violation
##----- DGP: Continuous time competing-risks with three encodings (CR, net-risk A/B)
##----- Strategies: Naive, Equi-confounding (EQC), Proximal Inference (PCI)
##----- Estimators per strategy: pooled logistic, Cox scalar, Cox time-varying
##----- Outputs: causal log HR (beta_2A) and cumulative incidence (ATT)
##----- FUNCTIONS (sourced by sim_unified_run.R)


library(data.table)
library(survival)
library(splines)


# Parameters and shared spline/CV basis ------------------------------------


default_params <- list(
  tau        = 53,                   # administrative follow-up (weeks)
  lam20      = 0.0019,               # test-positive baseline weekly rate  (~8% test >=1 positive)
  lam10      = 0.0032,               # test-negative baseline weekly rate  (~18% test >=1 negative)
  lambda_C0  = 0.004,                # censoring baseline rate
  b2A        = -0.7,                 # TRUE causal log hazard-ratio (constant across scenarios)
  b2X        = 0.3, b1X = 0.4,       # covariate effects (may differ across the two outcomes)
  b2U        = 0.8,                  # frailty coeff for positive outcome (set per scenario)
  b1U        = 0.8,                  # frailty coeff for negative outcome (set per scenario)
  b_CA       = 0.2, b_CX = 0.1,      # censoring depends on (A, X) only
  a0         = -0.2, aX = 0.5, aU = 1.0,  # treatment assignment (aU != 0 => U confounds A)
  gamma_A    = 0.5, gamma_Z = 0.6, gamma_X = 0.3, sigma_eps = 0.3  # U ~ gamma_A A + gamma_Z Z + gamma_X X + eps
)

## time basis shared by all pooled models
ns_knots  <- c(13, 26, 39)
ns_bknots <- c(1, 53)
## time-bin cut points for piecewise time-varying Cox
tv_cuts   <- c(13, 26, 39)


# Data-generating process -------------------------------------------------


gen_data <- function(N, p = default_params, scenario = "eqc_holds") {
  ## Scenario determines whether equi-confounding holds
  if (scenario == "eqc_violated") {
    p$b2U <- 1.0  # positive outcome more heavily confounded
    p$b1U <- 0.6
  } else {
    p$b2U <- 0.8
    p$b1U <- 0.8  # same as positive => equi-confounding holds
  }

  ## Covariates and unmeasured confounder
  X <- rnorm(N)
  Z <- rnorm(N)
  eps <- rnorm(N, 0, p$sigma_eps)
  U <- p$gamma_Z * Z + p$gamma_X * X + eps  # U will be re-computed conditional on A below
  
  ## Confounded treatment assignment
  A <- rbinom(N, 1, plogis(p$a0 + p$aX * X + p$aU * U))
  
  ## Recompute U conditional on A (reverse-order location shift for frailty bridge)
  U <- p$gamma_A * A + p$gamma_Z * Z + p$gamma_X * X + eps

  ## Cause-specific hazard rates
  rate2 <- p$lam20 * exp(p$b2A * A + p$b2X * X + p$b2U * U)  # test-positive (depends on A)
  rate1 <- p$lam10 * exp(              p$b1X * X + p$b1U * U)  # test-negative (NO A: valid NCO)
  crate  <- p$lambda_C0 * exp(p$b_CA * A + p$b_CX * X)       # LTFU (depends on A, X only)

  ## Generate event times
  C  <- pmin(p$tau, rexp(N, crate))
  T2 <- rexp(N, rate2)                                         # first positive test

  ## Recurrent negative tests: Poisson process on [0, C]
  nneg <- rpois(N, rate1 * C)
  negs <- data.table(
    id = rep(seq_len(N), nneg),
    t  = unlist(lapply(seq_len(N), function(i) if (nneg[i]) sort(runif(nneg[i], 0, C[i])) else numeric(0)))
  )

  ## First negative time
  T1f <- negs[, .(T1f = min(t)), by = id]
  dd  <- data.table(id = seq_len(N), X = X, Z = Z, U = U, A = A, C = C, T2 = T2, rate1 = rate1)
  dd  <- merge(dd, T1f, by = "id", all.x = TRUE)
  dd[is.na(T1f), T1f := Inf]

  list(dd = dd, negs = negs, p = p, scenario = scenario)
}


# Person-week expansion helpers (pooled discrete-time) --------------------


expand_terminal <- function(dd, end_wk, ev_wk) {
  reps <- pmax(as.integer(end_wk), 1L)
  out  <- dd[rep(seq_len(nrow(dd)), reps), .(id, X, Z, A)]
  out[, time := sequence(reps)]
  out[, .ev := rep(ev_wk, reps)]
  out[, event := as.integer(!is.na(.ev) & time == .ev)]
  out[, .ev := NULL]
  out[]
}

## Positive person-weeks on the net risk set {T2 > t}
pooled_pos_net <- function(dd, tau) {
  end_wk <- ceiling(pmin(dd$T2, dd$C, tau))
  ev_wk  <- ifelse(dd$T2 <= dd$C & dd$T2 <= tau, ceiling(dd$T2), NA_integer_)
  expand_terminal(dd, end_wk, ev_wk)
}

## Option B negatives: recurrent on the same risk set as positive {T2 > t}
pooled_neg_B <- function(dd, negs, tau) {
  base <- pooled_pos_net(dd, tau)
  base[, event := 0L]
  end  <- dd[, .(id, cap = pmin(T2, C, tau))]
  nw   <- merge(negs, end, by = "id")
  nw   <- nw[t < cap, .(id, time = ceiling(t))]
  nw   <- unique(nw)
  base[nw, on = .(id, time), event := 1L]
  base[]
}

## Option A negatives: first-negative time on its own risk set {T1 > t}
pooled_neg_A <- function(dd, tau) {
  end_wk <- ceiling(pmin(dd$T1f, dd$C, tau))
  ev_wk  <- ifelse(dd$T1f <= dd$C & dd$T1f <= tau, ceiling(dd$T1f), NA_integer_)
  expand_terminal(dd, end_wk, ev_wk)
}

## Competing-risk: first of {first negative, positive} is terminal on {min > t}
pooled_cr <- function(dd, tau) {
  crE     <- pmin(dd$T1f, dd$T2, dd$C, tau)
  end_wk  <- ceiling(crE)
  is_pos  <- dd$T2 <= dd$T1f & dd$T2 <= dd$C & dd$T2 <= tau
  is_neg  <- dd$T1f <  dd$T2 & dd$T1f <= dd$C & dd$T1f <= tau
  pos_ev  <- ifelse(is_pos, ceiling(dd$T2), NA_integer_)
  neg_ev  <- ifelse(is_neg, ceiling(dd$T1f), NA_integer_)
  list(pos = expand_terminal(dd, end_wk, pos_ev),
       neg = expand_terminal(dd, end_wk, neg_ev))
}


# Cox survival-object builders -----------------------------------------------


cox_pos_net <- function(dd, tau) {
  data.table(id = dd$id, X = dd$X, Z = dd$Z, A = dd$A,
             stop  = pmin(dd$T2, dd$C, tau),
             event = as.integer(dd$T2 <= dd$C & dd$T2 <= tau))
}

cox_neg_A <- function(dd, tau) {
  data.table(id = dd$id, X = dd$X, Z = dd$Z, A = dd$A,
             stop  = pmin(dd$T1f, dd$C, tau),
             event = as.integer(dd$T1f <= dd$C & dd$T1f <= tau))
}

## Andersen-Gill counting-process for recurrent negatives on {T2 > t}
cox_neg_B <- function(dd, negs, tau) {
  end <- dd[, .(id, X, Z, A, cap = pmin(T2, C, tau))]
  b   <- merge(negs, end, by = "id")
  b   <- b[t < cap][order(id, t)]
  b[, start := shift(t, fill = 0), by = id]
  ev  <- b[, .(id, X, Z, A, start, stop = t, event = 1L)]
  last <- b[, .(laststop = max(t)), by = id]
  fin <- merge(end, last, by = "id", all.x = TRUE)
  fin[is.na(laststop), laststop := 0]
  fin <- fin[cap > laststop, .(id, X, Z, A, start = laststop, stop = cap, event = 0L)]
  rbind(ev, fin)[order(id, start)]
}

cox_cr <- function(dd, tau) {
  crE    <- pmin(dd$T1f, dd$T2, dd$C, tau)
  is_pos <- dd$T2 <= dd$T1f & dd$T2 <= dd$C & dd$T2 <= tau
  is_neg <- dd$T1f <  dd$T2 & dd$T1f <= dd$C & dd$T1f <= tau
  list(
    pos = data.table(id = dd$id, X = dd$X, Z = dd$Z, A = dd$A, stop = crE, event = as.integer(is_pos)),
    neg = data.table(id = dd$id, X = dd$X, Z = dd$Z, A = dd$A, stop = crE, event = as.integer(is_neg))
  )
}


# Estimators: De-biased b2A by strategy (log-HR scale) ----------------------


## STRATEGY 1: NAIVE (no de-biasing, positive only)
est_naive_pooled <- function(long_pos, tau) {
  f <- event ~ ns(time, knots = ns_knots, Boundary.knots = ns_bknots) * A + X
  m2 <- suppressWarnings(glm(f, data = long_pos, family = binomial()))
  grid <- data.table(time = seq_len(tau), X = 0)
  d0 <- copy(grid)[, A := 0]; d1 <- copy(grid)[, A := 1]
  c2 <- predict(m2, d1, type = "link") - predict(m2, d0, type = "link")
  list(b2A = mean(c2), m1 = NULL, m2 = m2)  # m1 is NULL for naive (no NCO)
}

est_naive_cox_scalar <- function(surv_pos) {
  f <- coxph(Surv(stop, event) ~ A + X, data = surv_pos, ties = "efron")
  unname(coef(f)["A"])
}

est_naive_cox_tv <- function(surv_pos) {
  sp <- survSplit(Surv(stop, event) ~ ., data = as.data.frame(surv_pos),
                  cut = tv_cuts, episode = "bin", start = "tstart", end = "tstop")
  f  <- coxph(Surv(tstart, tstop, event) ~ A:factor(bin) + X, data = sp, ties = "efron")
  cf <- coef(f)[grep("A:factor", names(coef(f)))]
  mean(cf)
}


## STRATEGY 2: EQUI-CONFOUNDING (EQC) de-biasing = DiD on link scale
est_eqc_pooled <- function(long_pos, long_neg, tau) {
  f <- event ~ ns(time, knots = ns_knots, Boundary.knots = ns_bknots) * A + X
  m2 <- suppressWarnings(glm(f, data = long_pos, family = binomial()))
  m1 <- suppressWarnings(glm(f, data = long_neg, family = binomial()))
  grid <- data.table(time = seq_len(tau), X = 0)
  d0 <- copy(grid)[, A := 0]; d1 <- copy(grid)[, A := 1]
  c2 <- predict(m2, d1, type = "link") - predict(m2, d0, type = "link")
  c1 <- predict(m1, d1, type = "link") - predict(m1, d0, type = "link")
  list(b2A = mean(c2 - c1), m1 = m1, m2 = m2)
}

est_eqc_cox_scalar <- function(surv_pos, surv_neg, recurrent_neg) {
  f_pos <- coxph(Surv(stop, event) ~ A + X, data = surv_pos, ties = "efron")
  if (recurrent_neg) {
    f_neg <- coxph(Surv(start, stop, event) ~ A + X + cluster(id), data = surv_neg, ties = "efron")
  } else {
    f_neg <- coxph(Surv(stop, event) ~ A + X, data = surv_neg, ties = "efron")
  }
  unname(coef(f_pos)["A"] - coef(f_neg)["A"])
}

est_eqc_cox_tv <- function(surv_pos, surv_neg, recurrent_neg) {
  if (recurrent_neg) {
    sp_pos <- survSplit(Surv(stop, event) ~ ., data = as.data.frame(surv_pos),
                        cut = tv_cuts, episode = "bin", start = "tstart", end = "tstop")
    f_pos  <- coxph(Surv(tstart, tstop, event) ~ A:factor(bin) + X, data = sp_pos, ties = "efron")
    sp_neg <- survSplit(Surv(start, stop, event) ~ ., data = as.data.frame(surv_neg),
                        cut = tv_cuts, episode = "bin", start = "tstart", end = "tstop")
    f_neg  <- coxph(Surv(tstart, tstop, event) ~ A:factor(bin) + X + cluster(id), data = sp_neg, ties = "efron")
  } else {
    sp_pos <- survSplit(Surv(stop, event) ~ ., data = as.data.frame(surv_pos),
                        cut = tv_cuts, episode = "bin", start = "tstart", end = "tstop")
    f_pos  <- coxph(Surv(tstart, tstop, event) ~ A:factor(bin) + X, data = sp_pos, ties = "efron")
    sp_neg <- survSplit(Surv(stop, event) ~ ., data = as.data.frame(surv_neg),
                        cut = tv_cuts, episode = "bin", start = "tstart", end = "tstop")
    f_neg  <- coxph(Surv(tstart, tstop, event) ~ A:factor(bin) + X, data = sp_neg, ties = "efron")
  }
  
  cfp <- coef(f_pos)[grep("A:factor", names(coef(f_pos)))]
  cfn <- coef(f_neg)[grep("A:factor", names(coef(f_neg)))]
  common <- intersect(names(cfp), names(cfn))
  mean(cfp[common] - cfn[common])
}


## STRATEGY 3: PROXIMAL INFERENCE (PCI)
## Stage 1: Fit negative outcome on (A, Z, X) -> predict fitted link
## Stage 2: Fit positive outcome on (A, X, bridge) -> extract A coefficient
est_pci_pooled <- function(long_pos, long_neg, tau) {
  f_neg <- event ~ ns(time, knots = ns_knots, Boundary.knots = ns_bknots) * A + Z + X
  f_pos <- event ~ ns(time, knots = ns_knots, Boundary.knots = ns_bknots) * A + X
  
  m1 <- suppressWarnings(glm(f_neg, data = long_neg, family = binomial()))
  
  ## Create bridge covariate
  long_pos_aug <- copy(long_pos)
  long_pos_aug[, bridge := predict(m1, long_pos_aug, type = "link")]
  
  ## Fit positive outcome with bridge
  f_pos_bridge <- event ~ ns(time, knots = ns_knots, Boundary.knots = ns_bknots) * A + X + bridge
  m2 <- suppressWarnings(glm(f_pos_bridge, data = long_pos_aug, family = binomial()))
  
  grid <- data.table(time = seq_len(tau), X = 0, bridge = 0)
  d0 <- copy(grid)[, A := 0]; d1 <- copy(grid)[, A := 1]
  c2 <- predict(m2, d1, type = "link") - predict(m2, d0, type = "link")
  
  # Return both models (m1 is needed for CIF calculation)
  list(b2A = mean(c2), m1 = m1, m2 = m2, m1_for_cif = m1)
}

est_pci_cox_scalar <- function(surv_pos, surv_neg, recurrent_neg) {
  if (recurrent_neg) {
    f_neg <- coxph(Surv(start, stop, event) ~ A + Z + X + cluster(id), data = surv_neg, ties = "efron")
  } else {
    f_neg <- coxph(Surv(stop, event) ~ A + Z + X, data = surv_neg, ties = "efron")
  }
  
  ## Add bridge covariate (fitted log-hazard from stage 1)
  surv_pos_aug <- copy(surv_pos)
  if (recurrent_neg) {
    surv_pos_aug[, bridge := predict(f_neg, newdata = data.frame(A = A, Z = Z, X = X, start = 0, stop = stop, event = event), type = "lp")]
  } else {
    surv_pos_aug[, bridge := predict(f_neg, newdata = data.frame(A = A, Z = Z, X = X, stop = stop, event = event), type = "lp")]
  }
  
  f_pos <- coxph(Surv(stop, event) ~ A + X + bridge, data = surv_pos_aug, ties = "efron")
  unname(coef(f_pos)["A"])
}

est_pci_cox_tv <- function(surv_pos, surv_neg, recurrent_neg) {
  if (recurrent_neg) {
    f_neg <- coxph(Surv(start, stop, event) ~ A + Z + X + cluster(id), data = surv_neg, ties = "efron")
  } else {
    f_neg <- coxph(Surv(stop, event) ~ A + Z + X, data = surv_neg, ties = "efron")
  }
  
  surv_pos_aug <- copy(surv_pos)
  if (recurrent_neg) {
    surv_pos_aug[, bridge := predict(f_neg, newdata = data.frame(A = A, Z = Z, X = X, start = 0, stop = stop, event = event), type = "lp")]
  } else {
    surv_pos_aug[, bridge := predict(f_neg, newdata = data.frame(A = A, Z = Z, X = X, stop = stop, event = event), type = "lp")]
  }
  
  sp_pos <- survSplit(Surv(stop, event) ~ ., data = as.data.frame(surv_pos_aug),
                      cut = tv_cuts, episode = "bin", start = "tstart", end = "tstop")
  f_pos  <- coxph(Surv(tstart, tstop, event) ~ A:factor(bin) + X + bridge, data = sp_pos, ties = "efron")
  
  cf <- coef(f_pos)[grep("A:factor", names(coef(f_pos)))]
  mean(cf)
}


# CIF calculations -------------------------------------------------------


## Pooled CIF: net 1-KM or competing-risk AJ, de-biased at tau
pooled_cif <- function(m1, m2, dd, tau, competing, strategy = "eqc") {
  base <- CJ(id = dd$id, time = seq_len(tau))
  base <- merge(base, dd[, .(id, X, Z, A)], by = "id")
  
  f0 <- copy(base)[, A := 0]; f1 <- copy(base)[, A := 1]
  clamp <- function(x) pmin(pmax(x, 0), 1)
  
  if (strategy == "naive") {
    ## Naive: only positive model, A effect vs A=0
    l2_1 <- clamp(predict(m2, f1, type = "response")); l2_0 <- clamp(predict(m2, f0, type = "response"))
    D <- base[, .(id, time, A)]
    D[, `:=`(l2_1 = l2_1, l2_0 = l2_0)]
    if (!competing) {
      D[, `:=`(S1 = cumprod(1 - l2_1), S0 = cumprod(1 - l2_0)), by = id]
      D[time == tau & A == 1, .(risk1 = 1 - mean(S1), risk0 = 1 - mean(S0), risk0corr = 1 - mean(S0))]  # naive has no correction
    } else {
      # For competing risks, use only positive hazard (negative assumed constant at observed mean)
      D[, `:=`(Sall_1 = cumprod(pmax(1 - l2_1, 0)), Sall_0 = cumprod(pmax(1 - l2_0, 0))), by = id]
      D[, `:=`(Slag_1 = shift(Sall_1, fill = 1), Slag_0 = shift(Sall_0, fill = 1)), by = id]
      D[, `:=`(F1 = cumsum(l2_1 * Slag_1), F0 = cumsum(l2_0 * Slag_0)), by = id]
      D[time == tau & A == 1, .(risk1 = mean(F1), risk0 = mean(F0), risk0corr = mean(F0))]
    }
  } else if (strategy == "eqc") {
    ## EQC: DiD correction via negative model
    l2_0 <- clamp(predict(m2, f0, type = "response")); l2_1 <- clamp(predict(m2, f1, type = "response"))
    l1_0 <- clamp(predict(m1, f0, type = "response")); l1_1 <- clamp(predict(m1, f1, type = "response"))
    l2_c <- clamp(l2_0 * (l1_1 / l1_0))  # corrected untreated (multiplicative EQC)
    
    D <- base[, .(id, time, A)]
    D[, `:=`(l2_0 = l2_0, l2_1 = l2_1, l1_0 = l1_0, l1_1 = l1_1, l2_c = l2_c)]
    
    if (!competing) {
      D[, `:=`(S1 = cumprod(1 - l2_1), S0 = cumprod(1 - l2_0), Sc = cumprod(1 - l2_c)), by = id]
      D[time == tau & A == 1, .(risk1 = 1 - mean(S1), risk0 = 1 - mean(S0), risk0corr = 1 - mean(Sc))]
    } else {
      D[, `:=`(Sall_1 = cumprod(pmax(1 - l1_1 - l2_1, 0)),
               Sall_0 = cumprod(pmax(1 - l1_0 - l2_0, 0)),
               Sall_c = cumprod(pmax(1 - l1_1 - l2_c, 0))), by = id]
      D[, `:=`(Slag_1 = shift(Sall_1, fill = 1),
               Slag_0 = shift(Sall_0, fill = 1),
               Slag_c = shift(Sall_c, fill = 1)), by = id]
      D[, `:=`(F1 = cumsum(l2_1 * Slag_1),
               F0 = cumsum(l2_0 * Slag_0),
               Fc = cumsum(l2_c * Slag_c)), by = id]
      D[time == tau & A == 1, .(risk1 = mean(F1), risk0 = mean(F0), risk0corr = mean(Fc))]
    }
  } else if (strategy == "pci") {
    ## PCI: bridge is a proxy for the individual's fixed unmeasured confounder level and
    ## must be evaluated at their OBSERVED (A, Z, X) -- held fixed across the counterfactual
    ## contrast -- per the switching relation lambda_2^a = lambda_2(obs) * exp{beta_2A(a-A_obs)}.
    obs_bridge <- predict(m1, base, type = "link")
    f0$bridge <- obs_bridge; f1$bridge <- obs_bridge
    l2_1 <- clamp(predict(m2, f1, type = "response")); l2_0 <- clamp(predict(m2, f0, type = "response"))
    l1_0 <- clamp(predict(m1, f0, type = "response")); l1_1 <- clamp(predict(m1, f1, type = "response"))
    l2_c <- l2_0  # bridge-adjusted counterfactual risk already de-biased
    
    D <- base[, .(id, time, A)]
    D[, `:=`(l2_0 = l2_0, l2_1 = l2_1, l1_0 = l1_0, l1_1 = l1_1, l2_c = l2_c)]
    
    if (!competing) {
      D[, `:=`(S1 = cumprod(1 - l2_1), S0 = cumprod(1 - l2_0), Sc = cumprod(1 - l2_c)), by = id]
      D[time == tau & A == 1, .(risk1 = 1 - mean(S1), risk0 = 1 - mean(S0), risk0corr = 1 - mean(Sc))]
    } else {
      D[, `:=`(Sall_1 = cumprod(pmax(1 - l1_1 - l2_1, 0)),
               Sall_0 = cumprod(pmax(1 - l1_0 - l2_0, 0)),
               Sall_c = cumprod(pmax(1 - l1_1 - l2_c, 0))), by = id]
      D[, `:=`(Slag_1 = shift(Sall_1, fill = 1),
               Slag_0 = shift(Sall_0, fill = 1),
               Slag_c = shift(Sall_c, fill = 1)), by = id]
      D[, `:=`(F1 = cumsum(l2_1 * Slag_1),
               F0 = cumsum(l2_0 * Slag_0),
               Fc = cumsum(l2_c * Slag_c)), by = id]
      D[time == tau & A == 1, .(risk1 = mean(F1), risk0 = mean(F0), risk0corr = mean(Fc))]
    }
  }
}


## CIF CURVE version: returns risk1/risk0/risk0corr at EVERY week 1..tau (for quantile-band figure)
## Only implemented for the CR (competing-risks) encoding, pooled logistic estimator
pooled_cif_curve <- function(m1, m2, dd, tau, strategy = "eqc") {
  base <- CJ(id = dd$id, time = seq_len(tau))
  base <- merge(base, dd[, .(id, X, Z, A)], by = "id")
  f0 <- copy(base)[, A := 0]; f1 <- copy(base)[, A := 1]
  clamp <- function(x) pmin(pmax(x, 0), 1)

  if (strategy == "naive") {
    l2_1 <- clamp(predict(m2, f1, type = "response")); l2_0 <- clamp(predict(m2, f0, type = "response"))
    D <- base[, .(id, time, A)]
    D[, `:=`(l2_1 = l2_1, l2_0 = l2_0)]
    D[, `:=`(Sall_1 = cumprod(pmax(1 - l2_1, 0)), Sall_0 = cumprod(pmax(1 - l2_0, 0))), by = id]
    D[, `:=`(Slag_1 = shift(Sall_1, fill = 1), Slag_0 = shift(Sall_0, fill = 1)), by = id]
    D[, `:=`(F1 = cumsum(l2_1 * Slag_1), F0 = cumsum(l2_0 * Slag_0)), by = id]
    out <- D[A == 1, .(risk1 = mean(F1), risk0corr = mean(F0)), by = time]
  } else if (strategy == "eqc") {
    l2_0 <- clamp(predict(m2, f0, type = "response")); l2_1 <- clamp(predict(m2, f1, type = "response"))
    l1_0 <- clamp(predict(m1, f0, type = "response")); l1_1 <- clamp(predict(m1, f1, type = "response"))
    l2_c <- clamp(l2_0 * (l1_1 / l1_0))
    D <- base[, .(id, time, A)]
    D[, `:=`(l2_0 = l2_0, l2_1 = l2_1, l1_0 = l1_0, l1_1 = l1_1, l2_c = l2_c)]
    D[, `:=`(Sall_1 = cumprod(pmax(1 - l1_1 - l2_1, 0)),
             Sall_c = cumprod(pmax(1 - l1_1 - l2_c, 0))), by = id]
    D[, `:=`(Slag_1 = shift(Sall_1, fill = 1),
             Slag_c = shift(Sall_c, fill = 1)), by = id]
    D[, `:=`(F1 = cumsum(l2_1 * Slag_1), Fc = cumsum(l2_c * Slag_c)), by = id]
    out <- D[A == 1, .(risk1 = mean(F1), risk0corr = mean(Fc)), by = time]
  } else if (strategy == "pci") {
    ## bridge held fixed at each individual's OBSERVED (A, Z, X) -- see note in pooled_cif().
    obs_bridge <- predict(m1, base, type = "link")
    f0$bridge <- obs_bridge; f1$bridge <- obs_bridge
    l2_1 <- clamp(predict(m2, f1, type = "response")); l2_0 <- clamp(predict(m2, f0, type = "response"))
    l1_0 <- clamp(predict(m1, f0, type = "response")); l1_1 <- clamp(predict(m1, f1, type = "response"))
    D <- base[, .(id, time, A)]
    D[, `:=`(l2_0 = l2_0, l2_1 = l2_1, l1_0 = l1_0, l1_1 = l1_1)]
    D[, `:=`(Sall_1 = cumprod(pmax(1 - l1_1 - l2_1, 0)),
             Sall_c = cumprod(pmax(1 - l1_1 - l2_0, 0))), by = id]
    D[, `:=`(Slag_1 = shift(Sall_1, fill = 1),
             Slag_c = shift(Sall_c, fill = 1)), by = id]
    D[, `:=`(F1 = cumsum(l2_1 * Slag_1), Fc = cumsum(l2_0 * Slag_c)), by = id]
    out <- D[A == 1, .(risk1 = mean(F1), risk0corr = mean(Fc)), by = time]
  }
  out[]
}


## Oracle (true) untreated CIF curve for the CR encoding, standardized among the treated (ATT).
## Exact given the exponential/Poisson DGP: T2 ~ Exp(rate2), first-negative time ~ Exp(rate1),
## so the analytic competing-exponential-risks CIF applies at every week using each individual's
## TRUE (X, U) and the scenario's structural parameters (not estimated from a fitted model).
oracle_att_curve <- function(dd, p, tau) {
  trt  <- dd$A == 1
  r2_0 <- p$lam20 * exp(p$b2X * dd$X + p$b2U * dd$U)   # true untreated (A=0) rate2
  r1   <- dd$rate1                                       # true rate1 (already A-free)
  tt   <- seq_len(tau)
  F0 <- vapply(tt, function(t) {
    mean(((r2_0 / (r2_0 + r1)) * (1 - exp(-(r2_0 + r1) * t)))[trt])
  }, numeric(1))
  data.table(time = tt, F0_true = F0)
}


# One replication ---------------------------------------------------------


run_once <- function(N, p = default_params, scenario = "eqc_holds", seed = NULL, cox_encodings = "CR") {
  ## cox_encodings: which encodings to fit Cox (scalar/tv) estimators for (pooled is always fit for all).
  ## Default "CR" only, to save computation -- cross-estimator robustness is demonstrated on the
  ## primary (competing-risks) encoding; net-risk A/B use the pooled estimator only.
  if (!is.null(seed)) set.seed(seed)
  
  g  <- gen_data(N, p, scenario); dd <- g$dd; negs <- g$negs; tau <- p$tau
  
  ## Build pooled datasets for each encoding (always needed)
  P_pos  <- pooled_pos_net(dd, tau)
  P_negB <- pooled_neg_B(dd, negs, tau)
  P_negA <- pooled_neg_A(dd, tau)
  P_cr   <- pooled_cr(dd, tau)

  ## Results list
  res_b2A <- list()
  res_cif <- list()
  idx <- 1

  ## For each encoding (CR, netA, netB) x strategy (naive, eqc, pci) x estimator (pooled, [cox_scalar, cox_tv])
  
  for (encoding in c("CR", "netA", "netB")) {
    recurrent <- (encoding == "netB")  # netB uses recurrent negative events
    competing <- (encoding == "CR")
    fit_cox   <- encoding %in% cox_encodings
    
    ## Select the appropriate pooled datasets
    if (encoding == "CR") {
      Ppos <- P_cr$pos; Pneg <- P_cr$neg
    } else if (encoding == "netA") {
      Ppos <- P_pos; Pneg <- P_negA
    } else { # netB
      Ppos <- P_pos; Pneg <- P_negB
    }
    ## Build Cox survival objects only if needed for this encoding
    if (fit_cox) {
      if (encoding == "CR") {
        Xcr  <- cox_cr(dd, tau); Xpos <- Xcr$pos; Xneg <- Xcr$neg
      } else if (encoding == "netA") {
        Xpos <- cox_pos_net(dd, tau); Xneg <- cox_neg_A(dd, tau)
      } else {
        Xpos <- cox_pos_net(dd, tau); Xneg <- cox_neg_B(dd, negs, tau)
      }
    }

    ## NAIVE strategy (pooled only)
    naive_p <- est_naive_pooled(Ppos, tau)
    res_b2A[[idx]] <- data.table(encoding = encoding, strategy = "naive", estimator = "pooled",
                                  b2A = naive_p$b2A, b2A_true = p$b2A, scenario = scenario)
    res_cif[[length(res_cif) + 1]] <- cbind(encoding = encoding, strategy = "naive", estimator = "pooled",
                            pooled_cif(NULL, naive_p$m2, dd, tau, competing, "naive"),
                            scenario = scenario)
    idx <- idx + 1

    if (fit_cox) {
      res_b2A[[idx]] <- data.table(encoding = encoding, strategy = "naive", estimator = "cox_scalar",
                                    b2A = est_naive_cox_scalar(Xpos), b2A_true = p$b2A, scenario = scenario)
      idx <- idx + 1

      res_b2A[[idx]] <- data.table(encoding = encoding, strategy = "naive", estimator = "cox_tv",
                                    b2A = est_naive_cox_tv(Xpos), b2A_true = p$b2A, scenario = scenario)
      idx <- idx + 1
    }

    ## EQC strategy (pooled always; cox_scalar/cox_tv only for cox_encodings)
    eqc_p <- est_eqc_pooled(Ppos, Pneg, tau)
    res_b2A[[idx]] <- data.table(encoding = encoding, strategy = "eqc", estimator = "pooled",
                                  b2A = eqc_p$b2A, b2A_true = p$b2A, scenario = scenario)
    res_cif[[length(res_cif) + 1]] <- cbind(encoding = encoding, strategy = "eqc", estimator = "pooled",
                            pooled_cif(eqc_p$m1, eqc_p$m2, dd, tau, competing, "eqc"),
                            scenario = scenario)
    idx <- idx + 1

    if (fit_cox) {
      res_b2A[[idx]] <- data.table(encoding = encoding, strategy = "eqc", estimator = "cox_scalar",
                                    b2A = est_eqc_cox_scalar(Xpos, Xneg, recurrent), b2A_true = p$b2A, scenario = scenario)
      idx <- idx + 1

      res_b2A[[idx]] <- data.table(encoding = encoding, strategy = "eqc", estimator = "cox_tv",
                                    b2A = est_eqc_cox_tv(Xpos, Xneg, recurrent), b2A_true = p$b2A, scenario = scenario)
      idx <- idx + 1
    }

    ## PCI strategy (pooled always; cox_scalar/cox_tv only for cox_encodings)
    pci_p <- est_pci_pooled(Ppos, Pneg, tau)
    res_b2A[[idx]] <- data.table(encoding = encoding, strategy = "pci", estimator = "pooled",
                                  b2A = pci_p$b2A, b2A_true = p$b2A, scenario = scenario)
    res_cif[[length(res_cif) + 1]] <- cbind(encoding = encoding, strategy = "pci", estimator = "pooled",
                            pooled_cif(pci_p$m1, pci_p$m2, dd, tau, competing, "pci"),
                            scenario = scenario)
    idx <- idx + 1

    if (fit_cox) {
      res_b2A[[idx]] <- data.table(encoding = encoding, strategy = "pci", estimator = "cox_scalar",
                                    b2A = est_pci_cox_scalar(Xpos, Xneg, recurrent), b2A_true = p$b2A, scenario = scenario)
      idx <- idx + 1

      res_b2A[[idx]] <- data.table(encoding = encoding, strategy = "pci", estimator = "cox_tv",
                                    b2A = est_pci_cox_tv(Xpos, Xneg, recurrent), b2A_true = p$b2A, scenario = scenario)
      idx <- idx + 1
    }
  }

  list(b2A = rbindlist(res_b2A), cif = rbindlist(res_cif))
}
