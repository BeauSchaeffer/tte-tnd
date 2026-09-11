##----- Kaiser Causal TTE-TND
##----- Simulation study 2: EQC vs PCI under equi-confounding VIOLATION
##----- FUNCTIONS (sourced by sim_eqc_pci_run.R)
##-----
##----- Additive-hazard, Gaussian-frailty DGP so the proximal location shift is
##----- preserved under survivor conditioning => PCI is EXACTLY valid (additive case);
##----- a multiplicative-hazard variant (mirroring the main analysis) is also provided,
##----- where PCI is only APPROXIMATELY valid (survivor conditioning breaks the shift).
##----- Both include an NCE Z (excluded from the hazards) and toggle equi-confounding
##----- via beta_2U vs beta_1U:
##-----   scenario "eqc_holds"    : beta_1U = beta_2U  (EQC exact,  PCI valid)
##-----   scenario "eqc_violated" : beta_1U != beta_2U (EQC biased, PCI valid)
##----- The EQC difference returns beta_2A + gamma_A(beta_2U-beta_1U); PCI uses Z and is
##----- invariant to the equi-confounding gap.


library(data.table)
library(splines)


# Parameters --------------------------------------------------------------


pci_params <- function(scenario = c("eqc_holds", "eqc_violated"),
                       dgp = c("additive", "multiplicative")) {
  scenario <- match.arg(scenario); dgp <- match.arg(dgp)
  common <- list(tau = 53, a0 = -0.2, aX = 0.4, cens = 0.004, censA = 0.2, censX = 0.1,
                 scenario = scenario, dgp = dgp)
  if (dgp == "additive") {
    ## additive hazards + Gaussian frailty => location shift preserved => PCI EXACT
    b2U <- if (scenario == "eqc_holds") 0.015 else 0.027
    return(c(common, list(
      family = "gaussian", b2A = -0.02,
      gA = 0.4, gZ = 0.55, gX = 0.3, sU = 0.25,
      b10 = 0.06, b1X = 0.003, b1U = 0.015,
      b20 = 0.08, b2X = 0.002, b2U = b2U)))
  }
  ## multiplicative hazards (mirrors the main analysis) => PCI APPROXIMATE (survivor conditioning)
  b2U <- if (scenario == "eqc_holds") 0.6 else 1.0
  c(common, list(
    family = "binomial", b2A = -0.5,                    # true causal log hazard-ratio
    gA = 0.5, gZ = 0.6, gX = 0.3, sU = 0.3,
    b10 = 0.02, b1X = 0.3, b1U = 0.6,
    b20 = 0.006, b2X = 0.3, b2U = b2U))
}

ns_time <- function() "ns(time, knots = c(13, 26, 39), Boundary.knots = c(1, 53))"


# Data-generating process (additive hazards, competing risks) ------------


gen_pci <- function(N, p) {
  X <- rnorm(N); Z <- rnorm(N)
  A <- rbinom(N, 1, plogis(p$a0 + p$aX * X))
  U <- p$gA * A + p$gZ * Z + p$gX * X + rnorm(N, 0, p$sU)   # Gaussian location shift
  cl <- function(x) pmin(pmax(x, 0), 1)
  if (p$dgp == "additive") {
    lam1 <- cl(p$b10 + p$b1X * X + p$b1U * U)                       # test-negative (NCO)
    lam2 <- cl(p$b20 + p$b2A * A + p$b2X * X + p$b2U * U)           # test-positive
  } else {
    lam1 <- cl(p$b10 * exp(p$b1X * X + p$b1U * U))                  # multiplicative NCO
    lam2 <- cl(p$b20 * exp(p$b2A * A + p$b2X * X + p$b2U * U))      # multiplicative test-positive
  }

  Cw <- pmin(p$tau, ceiling(rexp(N, p$cens * exp(p$censA * A + p$censX * X))))
  atrisk <- rep(TRUE, N); Tt <- rep(NA_integer_, N); cause <- rep(0L, N)
  for (t in 1:p$tau) {
    idx <- which(atrisk); if (!length(idx)) break
    ev  <- runif(length(idx)) < (lam1[idx] + lam2[idx])
    hit <- idx[ev]
    pos <- runif(length(hit)) < (lam2[hit] / (lam1[hit] + lam2[hit]))
    Tt[hit] <- t
    cause[hit[pos]] <- 2L; cause[hit[!pos]] <- 1L
    atrisk[hit] <- FALSE
  }
  data.table(id = seq_len(N), X = X, Z = Z, A = A, U = U, Cw = Cw, Tt = Tt, cause = cause)
}

## person-week expansion, competing-risks encoding (first event terminal, risk set {T>t})
long_cr <- function(dd, tau) {
  evt    <- !is.na(dd$Tt) & dd$Tt <= dd$Cw               # event observed before censoring
  end    <- pmin(ifelse(is.na(dd$Tt), dd$Cw, pmin(dd$Tt, dd$Cw)), tau)
  ev_wk  <- ifelse(evt & dd$Tt <= tau, dd$Tt, NA_integer_)
  reps   <- pmax(as.integer(end), 1L)
  L <- dd[rep(seq_len(nrow(dd)), reps), .(id, X, Z, A)]
  L[, time := sequence(reps)]
  L[, ew := rep(ev_wk, reps)]
  L[, cz := rep(dd$cause, reps)]
  L[, Y_pos := as.integer(!is.na(ew) & time == ew & cz == 2L)]
  L[, Y_neg := as.integer(!is.na(ew) & time == ew & cz == 1L)]
  L[, c("ew", "cz") := NULL]
  L[]
}


# Estimators (additive / identity-link; A-contrast = hazard difference) --


## additive A-effect averaged over time (identity link => response-scale contrast)
Acon <- function(m, tau, extra = NULL) {
  g0 <- data.table(time = seq_len(tau), A = 0, X = 0)
  if (!is.null(extra)) g0[[extra]] <- 0
  g1 <- copy(g0)[, A := 1]
  mean(predict(m, g1) - predict(m, g0))
}

b2A_naive <- function(L, tau, fam) {
  f <- as.formula(sprintf("Y_pos ~ %s * A + X", ns_time()))
  Acon(glm(f, data = L, family = fam), tau)
}

## equi-confounding: test-positive minus test-negative treatment effect (link scale)
b2A_eqc <- function(L, tau, fam) {
  fp <- as.formula(sprintf("Y_pos ~ %s * A + X", ns_time()))
  fn <- as.formula(sprintf("Y_neg ~ %s * A + X", ns_time()))
  m2 <- glm(fp, data = L, family = fam)
  m1 <- glm(fn, data = L, family = fam)
  Acon(m2, tau) - Acon(m1, tau)
}

## proximal: stage-1 NCO hazard on (A,Z,X) is the bridge; stage-2 A effect = beta_2A
b2A_pci <- function(L, tau, fam) {
  f1 <- as.formula(sprintf("Y_neg ~ %s * (A + Z + X)", ns_time()))
  s1 <- glm(f1, data = L, family = fam)
  L <- copy(L)[, bridge := predict(s1, L)]              # link-scale bridge (log-hazard / hazard)
  f2 <- as.formula(sprintf("Y_pos ~ %s * A + X + bridge", ns_time()))
  s2 <- glm(f2, data = L, family = fam)
  Acon(s2, tau, extra = "bridge")
}


# One replication ---------------------------------------------------------


run_pci_once <- function(N, p) {
  fam <- if (p$family == "gaussian") gaussian() else binomial()
  dd <- gen_pci(N, p); L <- long_cr(dd, p$tau); tau <- p$tau
  data.table(
    dgp      = p$dgp,
    scenario = p$scenario,
    method   = c("naive", "eqc", "pci"),
    b2A      = c(b2A_naive(L, tau, fam), b2A_eqc(L, tau, fam), b2A_pci(L, tau, fam)),
    truth    = p$b2A
  )
}
