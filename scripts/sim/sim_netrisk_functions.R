##----- Kaiser Causal TTE-TND
##----- Simulation study: net-risk (Option A/B) vs competing-risk EQC de-biasing
##----- Estimators: pooled discrete-time logistic, Cox (time-varying coefs), Cox (scalar)
##----- FUNCTIONS (sourced by sim_netrisk_run.R)
##-----
##----- Data-generating process (continuous time, shared frailty):
##-----   U ~ N(0,1) shared frailty; A ~ Bern(expit(a0 + aX X + aU U))  (aU: confounding)
##-----   test-positive : first-passage hazard  rate2 = lam20 exp(b2A A + b2X X + bU U)
##-----   test-negative : recurrent Poisson rate rate1 = lam10 exp(       b1X X + bU U)   (NO A: NCO)
##-----   bU is SHARED by both hazards  => multiplicative equi-confounding holds exactly.
##----- The negatives are recurrent on [0, C]; each analysis FORMULATION uses them differently:
##-----   competing-risk : first event of {first negative, positive} is terminal (risk set {min>t})
##-----   net Option A   : positive on {T2>t}; NCO = first-negative time on its own set {T1A>t}
##-----   net Option B   : positive on {T2>t}; NCO = recurrent negative intensity on {T2>t}
##----- Truth: b2A is a known constant; frailty makes the MARGINAL hazards non-proportional,
##----- so scalar Cox (which collapses the time-varying marginal log-HR) is biased, while the
##----- pooled and time-varying-coefficient Cox estimators de-bias pointwise and are consistent.


library(data.table)
library(survival)
library(splines)


# Parameters --------------------------------------------------------------


default_params <- list(
  tau   = 53,                 # administrative follow-up (weeks)
  lam20 = 0.0019,             # test-positive baseline weekly rate  (~8% test >=1 positive)
  lam10 = 0.0032,             # test-negative baseline weekly rate  (~18% test >=1 negative)
  b2A   = -0.7,               # TRUE causal log hazard-ratio for testing positive (constant)
  b2X   = 0.3, b1X = 0.4,     # covariate effects (may differ across the two outcomes)
  bU    = 0.8,                # SHARED frailty coefficient  => equi-confounding
  a0 = -0.2, aX = 0.5, aU = 1.0,      # treatment assignment (aU != 0 => U confounds A)
  cens0 = 0.004, censA = 0.2, censX = 0.1  # LTFU rate depends on (A, X) only: C _||_ T | A, X
)

## time basis shared by the pooled models (fixed knots so predict() is stable)
ns_knots  <- c(13, 26, 39)
ns_bknots <- c(1, 53)
## time-bin cut points for the piecewise time-varying Cox
tv_cuts   <- c(13, 26, 39)


# Data-generating process -------------------------------------------------


gen_data <- function(N, p = default_params) {
  X <- rnorm(N); U <- rnorm(N)
  A <- rbinom(N, 1, plogis(p$a0 + p$aX * X + p$aU * U))

  rate2 <- p$lam20 * exp(p$b2A * A + p$b2X * X + p$bU * U)   # test-positive
  rate1 <- p$lam10 * exp(              p$b1X * X + p$bU * U)  # test-negative (no A)
  crate <- p$cens0 * exp(p$censA * A + p$censX * X)          # LTFU

  C  <- pmin(p$tau, rexp(N, crate))
  T2 <- rexp(N, rate2)

  ## recurrent negatives: Poisson(rate1) on [0, C]
  nneg <- rpois(N, rate1 * C)
  negs <- data.table(
    id = rep(seq_len(N), nneg),
    t  = unlist(lapply(seq_len(N), function(i) if (nneg[i]) sort(runif(nneg[i], 0, C[i])) else numeric(0)))
  )

  dd <- data.table(id = seq_len(N), X = X, U = U, A = A, C = C, T2 = T2, rate1 = rate1)
  T1 <- negs[, .(T1f = min(t)), by = id]
  dd <- merge(dd, T1, by = "id", all.x = TRUE)
  dd[is.na(T1f), T1f := Inf]
  list(dd = dd, negs = negs, p = p)
}


# Person-week expansion helpers (pooled discrete-time) --------------------


## expand a single terminal event to person-weeks 1..end_wk with event flag in ev_wk
expand_terminal <- function(dd, end_wk, ev_wk) {
  reps <- pmax(as.integer(end_wk), 1L)
  out  <- dd[rep(seq_len(nrow(dd)), reps), .(id, X, A)]
  out[, time := sequence(reps)]
  out[, .ev := rep(ev_wk, reps)]
  out[, event := as.integer(!is.na(.ev) & time == .ev)]
  out[, .ev := NULL]
  out[]
}

## positive person-weeks on the net risk set {T2 > t} (shared by Option A and B)
pooled_pos_net <- function(dd, tau) {
  end_wk <- ceiling(pmin(dd$T2, dd$C, tau))
  ev_wk  <- ifelse(dd$T2 <= dd$C & dd$T2 <= tau, ceiling(dd$T2), NA_integer_)
  expand_terminal(dd, end_wk, ev_wk)
}

## Option B negatives: recurrent markers on the SAME rows as pooled_pos_net
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

## Option A negatives: first-negative time on its own risk set {T1A > t}
pooled_neg_A <- function(dd, tau) {
  end_wk <- ceiling(pmin(dd$T1f, dd$C, tau))
  ev_wk  <- ifelse(dd$T1f <= dd$C & dd$T1f <= tau, ceiling(dd$T1f), NA_integer_)
  expand_terminal(dd, end_wk, ev_wk)
}

## competing-risk: first of {first negative, positive} is terminal on {min > t}
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


# Cox survival-object builders (continuous time) --------------------------


cox_pos_net <- function(dd, tau) {
  data.table(id = dd$id, X = dd$X, A = dd$A,
             stop  = pmin(dd$T2, dd$C, tau),
             event = as.integer(dd$T2 <= dd$C & dd$T2 <= tau))
}

cox_neg_A <- function(dd, tau) {
  data.table(id = dd$id, X = dd$X, A = dd$A,
             stop  = pmin(dd$T1f, dd$C, tau),
             event = as.integer(dd$T1f <= dd$C & dd$T1f <= tau))
}

## Andersen-Gill counting-process intervals for recurrent negatives on {T2 > t}
cox_neg_B <- function(dd, negs, tau) {
  end <- dd[, .(id, X, A, cap = pmin(T2, C, tau))]
  b   <- merge(negs, end, by = "id")
  b   <- b[t < cap][order(id, t)]
  b[, start := shift(t, fill = 0), by = id]
  ev  <- b[, .(id, X, A, start, stop = t, event = 1L)]
  last <- b[, .(laststop = max(t)), by = id]
  fin <- merge(end, last, by = "id", all.x = TRUE)
  fin[is.na(laststop), laststop := 0]
  fin <- fin[cap > laststop, .(id, X, A, start = laststop, stop = cap, event = 0L)]
  rbind(ev, fin)[order(id, start)]
}

cox_cr <- function(dd, tau) {
  crE    <- pmin(dd$T1f, dd$T2, dd$C, tau)
  is_pos <- dd$T2 <= dd$T1f & dd$T2 <= dd$C & dd$T2 <= tau
  is_neg <- dd$T1f <  dd$T2 & dd$T1f <= dd$C & dd$T1f <= tau
  list(
    pos = data.table(id = dd$id, X = dd$X, A = dd$A, stop = crE, event = as.integer(is_pos)),
    neg = data.table(id = dd$id, X = dd$X, A = dd$A, stop = crE, event = as.integer(is_neg))
  )
}


# Estimators: return de-biased b2A summary (log-HR scale) -----------------


## pooled discrete-time logistic; b2A(t) = pos A-contrast - neg A-contrast, averaged over t
est_pooled <- function(long_pos, long_neg, tau) {
  f <- event ~ ns(time, knots = ns_knots, Boundary.knots = ns_bknots) * A + X
  m2 <- suppressWarnings(glm(f, data = long_pos, family = binomial()))
  m1 <- suppressWarnings(glm(f, data = long_neg, family = binomial()))
  grid <- data.table(time = seq_len(tau), X = 0)
  d0 <- copy(grid)[, A := 0]; d1 <- copy(grid)[, A := 1]
  c2 <- predict(m2, d1, type = "link") - predict(m2, d0, type = "link")
  c1 <- predict(m1, d1, type = "link") - predict(m1, d0, type = "link")
  ## naive = positive A-contrast alone (no NCO de-biasing, U unadjusted)
  list(b2A = mean(c2 - c1), naive = mean(c2), b2A_t = c2 - c1, m1 = m1, m2 = m2)
}

## single-model A effect: scalar log-HR
cox_A_scalar <- function(surv, recurrent = FALSE) {
  f <- if (recurrent) coxph(Surv(start, stop, event) ~ A + X + cluster(id), data = surv, ties = "efron")
       else            coxph(Surv(stop, event) ~ A + X, data = surv, ties = "efron")
  unname(coef(f)["A"])
}

## single-model A effect: piecewise time-varying (named by time bin)
tv_binA <- function(fit) {
  cf <- coef(fit)
  cf <- cf[grep("A:factor\\(bin\\)", names(cf))]
  names(cf) <- sub(".*bin\\)", "", names(cf))   # keep bin index as name
  cf
}
cox_A_tv <- function(surv, recurrent = FALSE) {
  if (recurrent) {
    sp <- survSplit(Surv(start, stop, event) ~ ., data = as.data.frame(surv),
                    cut = tv_cuts, episode = "bin", start = "tstart", end = "tstop")
    f  <- coxph(Surv(tstart, tstop, event) ~ A:factor(bin) + X + cluster(id), data = sp, ties = "efron")
  } else {
    sp <- survSplit(Surv(stop, event) ~ ., data = as.data.frame(surv),
                    cut = tv_cuts, episode = "bin", start = "tstart", end = "tstop")
    f  <- coxph(Surv(tstart, tstop, event) ~ A:factor(bin) + X, data = sp, ties = "efron")
  }
  tv_binA(f)
}

## de-biased b2A = positive A-effect - negative (NCO) A-effect
est_cox_scalar <- function(surv_pos, surv_neg, recurrent_neg) {
  cox_A_scalar(surv_pos, FALSE) - cox_A_scalar(surv_neg, recurrent_neg)
}
est_cox_tv <- function(surv_pos, surv_neg, recurrent_neg) {
  bp <- cox_A_tv(surv_pos, FALSE); bn <- cox_A_tv(surv_neg, recurrent_neg)
  common <- intersect(names(bp), names(bn))
  mean(bp[common] - bn[common])
}


# Pooled CIF (secondary): net 1-KM (A/B) or competing-risk AJ (CR) --------


## per-person predicted hazards over weeks 1..L, standardized among the TREATED (ATT)
pooled_cif <- function(m1, m2, dd, tau, competing) {
  base <- CJ(id = dd$id, time = seq_len(tau))
  base <- merge(base, dd[, .(id, X, A)], by = "id")   # keep observed A for ATT standardization
  f0 <- copy(base)[, A := 0]; f1 <- copy(base)[, A := 1]
  clamp <- function(x) pmin(pmax(x, 0), 1)
  l2_0 <- clamp(predict(m2, f0, type = "response")); l2_1 <- clamp(predict(m2, f1, type = "response"))
  l1_0 <- clamp(predict(m1, f0, type = "response")); l1_1 <- clamp(predict(m1, f1, type = "response"))
  l2_c <- clamp(l2_0 * (l1_1 / l1_0))                 # corrected untreated (multiplicative EQC)

  D <- base[, .(id, time, A)]
  D[, `:=`(l2_0 = l2_0, l2_1 = l2_1, l1_0 = l1_0, l1_1 = l1_1, l2_c = l2_c)]
  if (!competing) {
    ## net: single-decrement 1 - prod(1 - lambda2); risk1 = observed treated (A=1), risk0corr = corrected A=0
    D[, `:=`(S1 = cumprod(1 - l2_1), S0 = cumprod(1 - l2_0), Sc = cumprod(1 - l2_c)), by = id]
    D[time == tau & A == 1, .(risk1 = 1 - mean(S1), risk0 = 1 - mean(S0), risk0corr = 1 - mean(Sc))]
  } else {
    ## competing-risk Aalen-Johansen: cumsum(lambda2 * S_lag), S all-cause (corrected uses lambda1 at A=1)
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


# Cox-based CIF (Breslow baseline -> per-week hazards -> same cumulation) --


## weekly baseline hazard increments (A=0, X=0) from a fitted Cox model
cox_baseline_weekly <- function(fit, tau) {
  bh <- basehaz(fit, centered = FALSE)
  H0 <- approx(c(0, bh$time), c(0, bh$hazard), xout = 0:tau, method = "constant", rule = 2)$y
  pmax(diff(H0), 0)
}

## per-week A log-HR (scalar: constant; tv: piecewise by time bin)
week_coefA <- function(fit, tv, tau) {
  if (!tv) return(rep(unname(coef(fit)["A"]), tau))
  binw <- findInterval(seq_len(tau) - 0.5, c(0, tv_cuts))   # 1..(len+1), matches survSplit episodes
  g <- tv_binA(fit)
  v <- as.numeric(g[as.character(binw)]); v[is.na(v)] <- 0
  v
}

## de-biased net (1-KM) / competing (AJ) CIF at week tau from Cox fits, among the TREATED (ATT)
cox_cif <- function(fit2, fit1, dd, tau, competing, tv) {
  N   <- nrow(dd); X <- dd$X; keep <- dd$A == 1
  h02 <- cox_baseline_weekly(fit2, tau)
  b2X <- unname(coef(fit2)["X"])
  b2A <- week_coefA(fit2, tv, tau)          # positive A effect (confounded)
  b1A <- week_coefA(fit1, tv, tau)          # NCO A effect (de-biasing term)
  emat <- function(v) matrix(v, N, tau, byrow = TRUE)

  base2  <- outer(exp(b2X * X), h02)         # lambda2(A=0, X): N x tau
  lam2_0 <- base2
  lam2_1 <- base2 * emat(exp(b2A))           # lambda2(A=1, X)
  lam2_c <- base2 * emat(exp(b1A))           # corrected untreated (multiplicative EQC)

  if (!competing) {
    f <- function(l) mean((1 - exp(-rowSums(l)))[keep])
    data.table(risk1 = f(lam2_1), risk0 = f(lam2_0), risk0corr = f(lam2_c))
  } else {
    h01 <- cox_baseline_weekly(fit1, tau); b1X <- unname(coef(fit1)["X"])
    base1  <- outer(exp(b1X * X), h01)
    lam1_0 <- base1; lam1_1 <- base1 * emat(exp(b1A))
    aj <- function(l2, l1) {
      tot  <- l1 + l2
      Slag <- exp(-(t(apply(tot, 1, cumsum)) - tot))   # survival to start of each week
      mean(rowSums(l2 * Slag)[keep])
    }
    data.table(risk1 = aj(lam2_1, lam1_1), risk0 = aj(lam2_0, lam1_0), risk0corr = aj(lam2_c, lam1_1))
  }
}

## fit the positive/negative Cox models (scalar or tv) then form the CIF
fit_cox_pos <- function(surv, tv) {
  if (!tv) return(coxph(Surv(stop, event) ~ A + X, data = surv, ties = "efron"))
  sp <- survSplit(Surv(stop, event) ~ ., data = as.data.frame(surv),
                  cut = tv_cuts, episode = "bin", start = "tstart", end = "tstop")
  coxph(Surv(tstart, tstop, event) ~ A:factor(bin) + X, data = sp, ties = "efron")
}
fit_cox_neg <- function(surv, recurrent, tv) {
  if (!tv) {
    if (recurrent) return(coxph(Surv(start, stop, event) ~ A + X + cluster(id), data = surv, ties = "efron"))
    return(coxph(Surv(stop, event) ~ A + X, data = surv, ties = "efron"))
  }
  if (recurrent) {
    sp <- survSplit(Surv(start, stop, event) ~ ., data = as.data.frame(surv),
                    cut = tv_cuts, episode = "bin", start = "tstart", end = "tstop")
    return(coxph(Surv(tstart, tstop, event) ~ A:factor(bin) + X + cluster(id), data = sp, ties = "efron"))
  }
  sp <- survSplit(Surv(stop, event) ~ ., data = as.data.frame(surv),
                  cut = tv_cuts, episode = "bin", start = "tstart", end = "tstop")
  coxph(Surv(tstart, tstop, event) ~ A:factor(bin) + X, data = sp, ties = "efron")
}
cif_cox <- function(surv_pos, surv_neg, recurrent, competing, dd, tau, tv) {
  cox_cif(fit_cox_pos(surv_pos, tv), fit_cox_neg(surv_neg, recurrent, tv), dd, tau, competing, tv)
}


# Oracle truth (closed form for exponential rates) -----------------------


oracle <- function(dd, p) {
  L <- p$tau; trt <- dd$A == 1                 # ATT: standardize among the treated
  r2 <- function(a) p$lam20 * exp(p$b2A * a + p$b2X * dd$X + p$bU * dd$U)
  r1 <- dd$rate1
  ## net first-positive CIF (single-decrement), standardized among the treated
  net <- function(a) mean((1 - exp(-r2(a) * L))[trt])
  ## crude competing-risk CIF with the first-negative (rate r1) competing
  crude <- function(a) mean(((r2(a) / (r2(a) + r1)) * (1 - exp(-(r2(a) + r1) * L)))[trt])
  list(b2A = p$b2A,
       net1 = net(1),   net0 = net(0),
       crude1 = crude(1), crude0 = crude(0))
}


# One replication ---------------------------------------------------------


run_once <- function(N, p = default_params) {
  g  <- gen_data(N, p); dd <- g$dd; negs <- g$negs; tau <- p$tau
  or <- oracle(dd, p)

  ## build datasets
  P_pos  <- pooled_pos_net(dd, tau)
  P_negB <- pooled_neg_B(dd, negs, tau)
  P_negA <- pooled_neg_A(dd, tau)
  P_cr   <- pooled_cr(dd, tau)

  X_pos  <- cox_pos_net(dd, tau)
  X_negB <- cox_neg_B(dd, negs, tau)
  X_negA <- cox_neg_A(dd, tau)
  X_cr   <- cox_cr(dd, tau)

  ## b2A recovery: 3 formulations x 3 estimators
  poolB <- est_pooled(P_pos,     P_negB,  tau)
  poolA <- est_pooled(P_pos,     P_negA,  tau)
  poolC <- est_pooled(P_cr$pos,  P_cr$neg, tau)

  res <- rbindlist(list(
    ## naive: net-T2 positive model only, adjusting X but NOT U and with NO NCO de-biasing
    data.table(formulation = "naive", estimator = "pooled",     b2A = poolB$naive),
    data.table(formulation = "naive", estimator = "cox_tv",     b2A = mean(cox_A_tv(X_pos, FALSE))),
    data.table(formulation = "naive", estimator = "cox_scalar", b2A = cox_A_scalar(X_pos, FALSE)),
    data.table(formulation = "netB", estimator = "pooled",     b2A = poolB$b2A),
    data.table(formulation = "netA", estimator = "pooled",     b2A = poolA$b2A),
    data.table(formulation = "CR",   estimator = "pooled",     b2A = poolC$b2A),
    data.table(formulation = "netB", estimator = "cox_tv",     b2A = est_cox_tv(X_pos, X_negB, TRUE)),
    data.table(formulation = "netA", estimator = "cox_tv",     b2A = est_cox_tv(X_pos, X_negA, FALSE)),
    data.table(formulation = "CR",   estimator = "cox_tv",     b2A = est_cox_tv(X_cr$pos, X_cr$neg, FALSE)),
    data.table(formulation = "netB", estimator = "cox_scalar", b2A = est_cox_scalar(X_pos, X_negB, TRUE)),
    data.table(formulation = "netA", estimator = "cox_scalar", b2A = est_cox_scalar(X_pos, X_negA, FALSE)),
    data.table(formulation = "CR",   estimator = "cox_scalar", b2A = est_cox_scalar(X_cr$pos, X_cr$neg, FALSE))
  ))
  res[, b2A_true := or$b2A]

  ## CIF by estimator: compare to formulation-specific oracle at landmark tau
  ## (cox_tv CIF omitted: basehaz has no baseline for the interaction-only TV model)
  mk <- function(form, est, cifdt, t1, t0)
    cbind(formulation = form, estimator = est, cifdt, truth1 = t1, truth0 = t0)
  cif <- rbindlist(list(
    mk("netB", "pooled",     pooled_cif(poolB$m1, poolB$m2, dd, tau, FALSE),           or$net1,   or$net0),
    mk("netA", "pooled",     pooled_cif(poolA$m1, poolA$m2, dd, tau, FALSE),           or$net1,   or$net0),
    mk("CR",   "pooled",     pooled_cif(poolC$m1, poolC$m2, dd, tau, TRUE),            or$crude1, or$crude0),
    mk("netB", "cox_scalar", cif_cox(X_pos,    X_negB,  TRUE,  FALSE, dd, tau, FALSE), or$net1,   or$net0),
    mk("netA", "cox_scalar", cif_cox(X_pos,    X_negA,  FALSE, FALSE, dd, tau, FALSE), or$net1,   or$net0),
    mk("CR",   "cox_scalar", cif_cox(X_cr$pos, X_cr$neg, FALSE, TRUE,  dd, tau, FALSE), or$crude1, or$crude0)
  ))
  list(b2A = res, cif = cif)
}
