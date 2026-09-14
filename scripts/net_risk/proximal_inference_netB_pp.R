##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Proximal Inference Analysis -- NET RISK, Option B
##----- Per-protocol, no censoring weights
##----- last updated 2026-09-14


# Packages ----------------------------------------------------------------


library(tidyverse)
library(tidycmprsk)
library(survival)
library(broom)
library(data.table)


# Data --------------------------------------------------------------------


### Net risk replaces the single three-level Y3 dataset with two panels
### that share a cohort but NOT a risk set:
###   Panel P (positive)  = data_Y2, follow-up to min(T2, C, tau)
###   Panel N (negative)  = recurrent negatives on {T2 > t}, built below
### Both are read from data_weekmatch.4/.

data_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch.4/"

data_Y2  <- read_rds(paste0(data_path, "data_Y2_weekmatch.rds"))
neg_hist <- read_rds(paste0(data_path, "neg_hist_weekmatch.rds"))

data_Y2 <- data_Y2 |>
  mutate(subclass = as.character(subclass))

### Input checks. Stage 1 and stage 2 are fit on different panels, so the
### bridge is only well defined if both cover exactly the same individuals.

mod_vars <- c("treatment",
              "sex_admin", "age_years", "bmi", "race", "charlson_cat_fac",
              "ndi", "prior_inf", "tests_count", "service_region",
              "last_vax_infect_weeks",
              "flu_vax")

stopifnot(!anyDuplicated(data_Y2$fake_mrn))
stopifnot(!anyNA(data_Y2 |> select(all_of(mod_vars))))
stopifnot(!anyNA(data_Y2$Y2_pp_t_trunc), !anyNA(data_Y2$Y2_pp_trunc))
stopifnot(all(data_Y2$Y2_pp_t_trunc > 0))
stopifnot(all(neg_hist$neg_t > 0))
stopifnot(all(neg_hist$fake_mrn %in% data_Y2$fake_mrn))   # no orphan negatives

message("Cohort: ", nrow(data_Y2), " individuals, ",
        n_distinct(data_Y2$subclass), " matched pairs")
message("Positive events (PP, truncated): ", sum(data_Y2$Y2_pp_trunc))

res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results_pp.5/"
dir.create(res_path, recursive = TRUE, showWarnings = FALSE)
stopifnot(dir.exists(res_path))


# Option B: Andersen-Gill expansion ---------------------------------------


### Identical construction to equi_confounding_netB_pp.R -- panel N is the
### recurrent negative-testing intensity on the positive panel's risk set.

### covariates carried onto every interval
ag_vars <- c("fake_mrn", "subclass", "treatment",
             # demog
             "sex_admin", "age_years", "bmi", "race", "charlson_cat_fac",
             # other
             "ndi", "prior_inf", "tests_count", "service_region",
             "last_vax_infect_weeks",
             # NEC
             "flu_vax")

ag_base <- data_Y2 |>
  select(all_of(ag_vars), cap = Y2_pp_t_trunc)

stopifnot(all(ag_base$cap > 0))

### negatives inside the PP risk set, one row per person-week
ag_neg <- neg_hist |>
  inner_join(ag_base |> select(fake_mrn, cap), by = "fake_mrn") |>
  filter(neg_t <= cap) |>                   # <= not <; whole weeks, so ties at the cap count
  distinct(fake_mrn, neg_t) |>              # collapse same-week repeats
  arrange(fake_mrn, neg_t)

### event intervals: (previous stop, neg_t] with event = 1
ag_event <- ag_neg |>
  group_by(fake_mrn) |>
  mutate(start = lag(neg_t, default = 0),
         stop  = neg_t,
         event = 1L) |>
  ungroup() |>
  select(fake_mrn, start, stop, event)

### trailing censored interval: (last neg_t, cap] with event = 0
### also covers everyone with no negatives at all, whose last stop is 0
ag_tail <- ag_neg |>
  group_by(fake_mrn) |>
  summarise(last_stop = max(neg_t), .groups = "drop") |>
  right_join(ag_base |> select(fake_mrn, cap), by = "fake_mrn") |>
  mutate(last_stop = coalesce(last_stop, 0)) |>
  filter(last_stop < cap) |>
  transmute(fake_mrn, start = last_stop, stop = cap, event = 0L)

ag_pp <- bind_rows(ag_event, ag_tail) |>
  inner_join(ag_base |> select(-cap), by = "fake_mrn") |>
  arrange(fake_mrn, start)

stopifnot(all(ag_pp$stop > ag_pp$start))                     # no zero-length intervals
stopifnot(n_distinct(ag_pp$fake_mrn) == nrow(data_Y2))       # nobody dropped
stopifnot(!anyNA(ag_pp |> select(all_of(mod_vars))))         # same completeness as panel P

### Person-time reconciliation -- intervals must tile [0, cap] exactly
ag_span <- ag_pp |>
  group_by(fake_mrn) |>
  summarise(span = sum(stop - start), .groups = "drop") |>
  inner_join(data_Y2 |> select(fake_mrn, cap = Y2_pp_t_trunc), by = "fake_mrn")

stopifnot(nrow(ag_span) == nrow(data_Y2))
stopifnot(isTRUE(all.equal(ag_span$span, ag_span$cap)))

message("AG expansion: ", nrow(ag_pp), " intervals, ",
        sum(ag_pp$event), " negative events in the {T2 > t} risk set")

setDT(ag_pp)
setDT(data_Y2)


# Cox model ---------------------------------------------------------------

# Per-protocol (no censoring weights)

### The bridge is the stage 1 linear predictor. Covariates in ag_pp are
### baseline-fixed, so the linear predictor is constant within person and
### maps cleanly from the interval-level stage 1 fit onto the person-level
### stage 2 panel. predict() needs the Surv() variables present in newdata
### even though type = "lp" does not use them, so they are supplied below.

bridge_newdata <- function(d) {
  d |>
    mutate(start = 0,
           stop  = Y2_pp_t_trunc,
           event = Y2_pp_trunc)
}

# Stage 1

prox_pp_s1 <- coxph(
  Surv(start, stop, event) ~ treatment +
    # demog
    sex_admin + age_years + bmi + race + charlson_cat_fac +
    # other
    ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
    # NEC
    flu_vax +
    cluster(subclass),
  data = ag_pp
)

# Predictions

data_Y2$p <- predict(prox_pp_s1, newdata = bridge_newdata(data_Y2), type = "lp")

stopifnot(!anyNA(data_Y2$p))

# Stage 2

prox_pp_s2 <- coxph(
  Surv(Y2_pp_t_trunc, Y2_pp_trunc) ~ treatment +
    # demog
    sex_admin + age_years + bmi + race + charlson_cat_fac +
    # other
    ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
    # predictions from S1
    p +
    # NO NEC
    cluster(subclass),
  data = data_Y2
)

stopifnot(prox_pp_s1$n == nrow(ag_pp))
stopifnot(prox_pp_s2$n == nrow(data_Y2))

netB.pci.pp.cox.pointest <- tidy(prox_pp_s2, conf.int = TRUE, exponentiate = TRUE)
saveRDS(netB.pci.pp.cox.pointest, paste0(res_path,"netB.pci.pp.cox.pointest.rds"))


# Bootstrap CIs for PP ---------------------------------------------------

num.boot <- 200

set.seed(1155)
seed <- floor(runif(num.boot)*10^8)

setkey(data_Y2, subclass)
setkey(ag_pp, subclass)

subclasses <- data_Y2[, unique(subclass)]
n_sub <- length(subclasses)

boot.results <- lapply(1:num.boot, function(i){
  
  t0 <- Sys.time()
  
  set.seed(seed[i])
  
  message("Starting PCI net-B PP Cox bootstrap ", i, " (seed=", seed[i], ")")
  
  # select matched pairs
  samp_sub <- sample(subclasses, size = n_sub, replace = TRUE)
  
  # build boot dataset efficiently via one join:
  # map draw index j -> sampled subclass, then join to replicate all rows per subclass
  # the SAME draw is applied to both panels, so stage 1 and stage 2 come from
  # one resample of pairs rather than two
  map <- data.table(j = seq_along(samp_sub), subclass = samp_sub)
  
  dat.boot <- data_Y2[map, on = "subclass", allow.cartesian = TRUE]
  ag.boot  <- ag_pp[map,   on = "subclass", allow.cartesian = TRUE]
  
  # new cluster id per draw
  dat.boot[, subclass_boot := paste0(subclass, ".", j)]
  ag.boot[,  subclass_boot := paste0(subclass, ".", j)]
  
  # run stage 1 model
  prox_pp_s1 <- coxph(
    Surv(start, stop, event) ~ treatment +
      # demog
      sex_admin + age_years + bmi + race + charlson_cat_fac +
      # other
      ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
      # NEC
      flu_vax +
      cluster(subclass_boot),
    data = ag.boot
  )
  
  # generate predictions
  
  dat.boot$p <- predict(prox_pp_s1, newdata = bridge_newdata(dat.boot), type = "lp")
  
  # run stage 2 model
  
  prox_pp_s2 <- coxph(
    Surv(Y2_pp_t_trunc, Y2_pp_trunc) ~ treatment +
      # demog
      sex_admin + age_years + bmi + race + charlson_cat_fac +
      # other
      ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
      # predictions from S1
      p +
      # NO NEC
      cluster(subclass_boot),
    data = dat.boot
  )
  
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  message("Finished bootstrap ", i, " in ", round(elapsed, 2), " minutes")
  
  return(c(
    sim=i,
    treatHR=unname(exp(prox_pp_s2$coefficients["treatment"]))
  ))
  
})

boot.long <- bind_rows(lapply(boot.results, tibble::as_tibble_row))
saveRDS(boot.long, paste0(res_path, "netB.pci.pp.cox.boot.long.rds"))
