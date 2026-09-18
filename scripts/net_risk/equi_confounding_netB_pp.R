##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Equi Confounding Analysis -- NET RISK, Option B
##----- Per-protocol, no censoring weights
##----- last updated 2026-09-14

# Packages ----------------------------------------------------------------


library(tidyverse)
library(tidycmprsk)
library(survival)
library(ggsurvfit)
library(riskRegression)
library(geepack)
library(data.table)


# Data --------------------------------------------------------------------


### Net risk replaces the single three-level Y3 dataset with two panels
### that share a cohort but NOT a risk set:
###   Panel P (positive)  = data_Y2, follow-up to min(T2, C, tau)
###   Panel N (negative)  = recurrent negatives on {T2 > t}, built below
### Both are read from data_weekmatch.3/.

data_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch.3/"

data_Y2  <- read_rds(paste0(data_path, "data_Y2_weekmatch.rds"))
neg_hist <- read_rds(paste0(data_path, "neg_hist_weekmatch.rds"))

data_Y2 <- data_Y2 |>
  mutate(subclass = as.character(subclass))

### Input checks. The two EQC models must be fit on the same individuals --
### if coxph silently drops rows for missingness in one panel but not the
### other, the hazard ratio contrasts two different cohorts.

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
message("Negative tests in neg_hist: ", nrow(neg_hist), " across ",
        n_distinct(neg_hist$fake_mrn), " individuals (",
        round(100 * n_distinct(neg_hist$fake_mrn) / nrow(data_Y2), 1), "% of cohort)")

res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results_pp.5/"
dir.create(res_path, recursive = TRUE, showWarnings = FALSE)
stopifnot(dir.exists(res_path))


# Option B: Andersen-Gill expansion ---------------------------------------


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

### Person-time reconciliation. The intervals must tile [0, cap] exactly, with
### no gaps or overlaps, so that the negative intensity is fit on precisely the
### same follow-up as the positive hazard. This is the property Option B rests
### on -- if it fails, the risk sets differ and every downstream number is wrong.

ag_span <- ag_pp |>
  group_by(fake_mrn) |>
  summarise(span = sum(stop - start), .groups = "drop") |>
  inner_join(data_Y2 |> select(fake_mrn, cap = Y2_pp_t_trunc), by = "fake_mrn")

stopifnot(nrow(ag_span) == nrow(data_Y2))
stopifnot(isTRUE(all.equal(ag_span$span, ag_span$cap)))

### Intervals must also start at 0 and be contiguous within person
ag_gap <- ag_pp |>
  group_by(fake_mrn) |>
  summarise(first_start = min(start),
            contiguous  = all(start[-1] == head(stop, -1)),
            .groups = "drop")

stopifnot(all(ag_gap$first_start == 0))
stopifnot(all(ag_gap$contiguous))

message("AG expansion: ", nrow(ag_pp), " intervals, ",
        sum(ag_pp$event), " negative events in the {T2 > t} risk set")

setDT(ag_pp)

# Sanity checks:
# ag_pp[, .(n_events = sum(event)), by = fake_mrn][, table(n_events)]


# Cox Model ---------------------------------------------------------------


# Per-protocol, no censoring weights

### fit1 is now the recurrent negative intensity (Andersen-Gill) rather
### than a first-test cause-specific hazard. The covariate set, the
### flu_vax NCE term, and the subclass clustering are unchanged.

eqc_pp_fit1 <- coxph(
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


eqc_pp_fit2 <- coxph(
  Surv(Y2_pp_t_trunc, Y2_pp_trunc) ~ treatment +
    # demog
    sex_admin + age_years + bmi + race + charlson_cat_fac +
    # other
    ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
    # NEC
    flu_vax +
    cluster(subclass),
  data = data_Y2
)

### Both fits must use the same individuals for the ratio to be meaningful.
### coxph drops incomplete cases silently; n in fit1 counts intervals, so
### compare distinct individuals rather than model n.

stopifnot(eqc_pp_fit2$n == nrow(data_Y2))
stopifnot(eqc_pp_fit1$n == nrow(ag_pp))

netB.eqc.pp.cox.pointest <- c(
  treatHR= unname( exp(eqc_pp_fit2$coefficients["treatment"]) / exp(eqc_pp_fit1$coefficients["treatment"]) ),
  fluvaxHR=unname( exp(eqc_pp_fit2$coefficients["flu_vax"]) / exp(eqc_pp_fit1$coefficients["flu_vax"]) )
)

saveRDS(netB.eqc.pp.cox.pointest, paste0(res_path,"netB.eqc.pp.cox.pointest.rds"))


# Bootstrap CIs for PP ---------------------------------------------------


num.boot <- 200

set.seed(1155)
seed <- floor(runif(num.boot)*10^8)

setDT(data_Y2)
setkey(data_Y2, subclass)
setkey(ag_pp, subclass)

subclasses <- data_Y2[, unique(subclass)]
n_sub <- length(subclasses)


boot.results <- lapply(1:num.boot, function(i){
  
  t0 <- Sys.time()
  
  set.seed(seed[i])
  
  message("Starting EQC net-B PP Cox bootstrap ", i, " (seed=", seed[i], ")")
  
  # select matched pairs
  samp_sub <- sample(subclasses, size = n_sub, replace = TRUE)
  
  # build boot dataset efficiently via one join:
  # map draw index j -> sampled subclass, then join to replicate all rows per subclass
  # the SAME draw is applied to both panels, so the two coefficients in
  # the ratio come from one resample of pairs rather than two
  map <- data.table(j = seq_along(samp_sub), subclass = samp_sub)
  
  dat.boot    <- data_Y2[map, on = "subclass", allow.cartesian = TRUE]
  ag.boot     <- ag_pp[map,   on = "subclass", allow.cartesian = TRUE]
  
  # new cluster id per draw
  dat.boot[, subclass_boot := paste0(subclass, ".", j)]
  ag.boot[,  subclass_boot := paste0(subclass, ".", j)]
  
  # run test neg model
  eqc_pp_fit1 <- coxph(
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
  
  # run test pos model
  eqc_pp_fit2 <- coxph(
    Surv(Y2_pp_t_trunc, Y2_pp_trunc) ~ treatment +
      # demog
      sex_admin + age_years + bmi + race + charlson_cat_fac +
      # other
      ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
      # NEC
      flu_vax +
      cluster(subclass_boot),
    data = dat.boot
  )
  
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  message("Finished bootstrap ", i, " in ", round(elapsed, 2), " minutes")
  
  return(c(
    sim=i,
    treatHR=unname(exp(eqc_pp_fit2$coefficients["treatment"])) / unname(exp(eqc_pp_fit1$coefficients["treatment"])),
    fluvaxHR=unname(exp(eqc_pp_fit2$coefficients["flu_vax"])) / unname(exp(eqc_pp_fit1$coefficients["flu_vax"]))
  ))
  
})

boot.long <- bind_rows(lapply(boot.results, tibble::as_tibble_row))
saveRDS(boot.long, paste0(res_path, "netB.eqc.pp.cox.boot.long.rds"))
