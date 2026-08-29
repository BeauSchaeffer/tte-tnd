##----- Beau Schaeffer
##----- Kaiser TTE-TND
##----- Data cleaning - recleaning to match on week rather than within enrollment
##----- includes edits for net risk estimation

# Packages ----------------------------------------------------------------

library(tidyverse)
library(MatchIt)
library(data.table)

# Load data ---------------------------------------------------------------

demog <- read_csv("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/std_elig_pop_250701.csv")
vax <- read_csv("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/vax_250701.csv")
labs <- read_csv("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/labs_250701.csv")

# Demographic / Patient Information ---------------------------------------

# demog |> select(fake_mrn) |> unique() |> nrow() # 950848

demog_clean <- demog |>
  
  filter(
    
    bmi<=100,
    # >100 (18 values)
    
    enroll_doses<=1,
    # 2 (1636 values), 3 (57), 4 (5), 5 (3)
    
    tests_count<=100,
    # >100 (62 values)
    
    charlson_cat != "No v", # 54133 No v
    # considering this as missing data
    
    last_enroll_booster=="NONE" | last_enroll_booster=="Primary series + 2 additional",
    # keep individuals who were P+1 or P+2 at the time of their last enrollment booster
    # I believe "NONE" indicates P+1
    # 894937
    
    months_last_vax_infect >= 3.0,
    #  last immunological event at least 3 months prior to enroll
    # redundant with last_pe_booster_weeks?
    # 871285
    
    last_pe_booster_weeks < -12.0
    # last booster at least 3 months prior to enroll
    # redundant with months_last_vax_infect?
    # 871285 after charlson filter
    
  ) |>
  
  mutate(
    
    sex_admin=factor(sex_admin),
    
    charlson_cat_fac=factor(charlson_cat,levels=c("0", "1", "2", "3+")), # see note above
    
    last_pe_booster=factor(last_pe_booster), # already restricted to 1 level
    
    last_enroll_booster=factor(last_enroll_booster, levels=c("NONE", "Primary series + 2 additional")), # see note above
    
    
    race=factor(race,levels=c("White", "Asian","Hispanic","Black","Unknown/ot")),
    
    
    service_region=factor(service_region,levels=c("Central valley", "East bay", "North bay",
                                                  "North valley", "Peninsula", "South bay")),
    
    
    charlson_cat_fac=factor(charlson_cat,levels=c("0", "1", "2", "3+"))) |>
  
  select(-enroll_infect, -enroll_infect_weeks)
# no longer considering infection analysis


# demog |> select(last_vax_infect_weeks) |> filter(last_vax_infect_weeks >= -12) |> nrow()
### 21184 with imm event within 12wk of 09/01
# demog_clean |> select(last_vax_infect_weeks) |> filter(last_vax_infect_weeks >= -12)

# demog |> select(months_last_vax_infect) |> summary()
# demog |> select(months_last_vax_infect) |> filter(months_last_vax_infect < 3) |> nrow()
### 25209 with imm event within 3mo of 09/01
### not relative dating weeks?

# demog |> select(last_pe_booster_weeks) |> summary()
# demog |> select(last_pe_booster_weeks) |> filter(last_pe_booster_weeks >= -12) |> nrow()
### 21184 with booster within 12wk of 09/01
# demog_clean |> select(last_pe_booster_weeks) |> summary()

# demog_clean |> count(fake_mrn) |> filter(n>1)
### no duplicate MRNs

# demog_clean |> select(fake_mrn) |> unique() |> nrow() # 871285
# 871285 --> this number should remain moving forward


### truncate follow up later in script


# Labs / Testing ----------------------------------------------------------


labs_clean <- labs |>
  
  filter(test_type=="SARS_COV_2_NAAT",
         pt_loc=="O",
         result!="no result") |>
  
  # mutate(BA = case_when(lab_weeks<0 ~ "before",
  #                       lab_weeks>= 0 ~ "after")) |>
  # select(BA) |> table()
  
  # before 1509552
  # after 265510
  
  filter(lab_weeks >= 0) # tests eligible as outcomes need positive value


# Join Demog + Labs -------------------------------------------------------


### For THREE level outcome (Pos, Neg, Cens)
### labs_clean_3L
### only keep the first lab result for individuals with multiple lab records
### positive or negative test considered an event

labs_clean_3L <- labs_clean |>
  select(fake_mrn, result, lab_weeks) |>
  group_by(fake_mrn) |>
  arrange(fake_mrn, lab_weeks) |> # added fake_mrn for visual check
  slice_min(lab_weeks, with_ties = F) |>
  ungroup() |>
  rename("result_3L" = "result",
         "lab_weeks_3L" = "lab_weeks")

# table(labs_clean_3L$result_3L)

### For TWO level outcome (Pos, Cens)
### labs_clean_2L
### only keep the first POS lab result for individuals with multiple lab records
### only positive test considered an event

labs_clean_2L <- labs_clean |>
  select(fake_mrn, result, lab_weeks) |>
  filter(result=="positive") |>
  group_by(fake_mrn) |>
  arrange(fake_mrn, lab_weeks) |> # added fake_mrn for visual check
  slice_min(lab_weeks, with_ties = F) |>
  ungroup() |>
  rename("result_2L" = "result",
         "lab_weeks_2L" = "lab_weeks")

# table(labs_clean_2L$result_2L)

#netrisk#
# retain all negative tests for net risk
labs_neg_all <- labs_clean |>
  filter(result == "negative") |>
  select(fake_mrn, lab_weeks) |>
  rename(neg_week = lab_weeks)
# labs_neg_all |> count(fake_mrn) |> pull(n) |> summary()
# labs_neg_all |> nrow()

demog_labs <- demog_clean |>
  left_join(labs_clean_3L, by = "fake_mrn") |>
  left_join(labs_clean_2L, by = "fake_mrn")

# demog_labs |> select(fake_mrn) |> unique() |> nrow() # 871285


# Vaccination -------------------------------------------------------------

# vax$vax_series |> table()

### 950848 unique fake_mrn

### 950848 Partial
### 950848 Complete primary series
### 950848 "Primary series + 1 additional"
### 414807 "P+2"
### 172609 "P+3"
### 5429 "P+4"
### 55 "P+5"
### 3235 "Other series"

vax_clean <- vax |>
  mutate(vax_name = factor(vax_name, levels=c("pfizer","moderna"))) |>
  
  filter(vax_series != "Other series") |>  # contain errors
  
  mutate(vax_series = factor(vax_series, levels=c("Partial", "Complete primary series",
                                                  "Primary series + 1 additional", "Primary series + 2 additional",
                                                  "Primary series + 3 additional", "Primary series + 4 additional",
                                                  "Primary series + 5 additional")))

### 414807 "P+2"
### check to see if bivalent
### 61130 not / 353677 bivalent

# vax_clean |> select(fake_mrn) |> unique() |> nrow() # 950848


# Join Demog + Labs + Vax -------------------------------------------------


demog_labs_vax <- demog_labs |>
  left_join(vax_clean, by="fake_mrn")

# demog_labs_vax |> select(fake_mrn) |> unique() |> nrow() # 871285

### creates multiple rows (per vax) per individual
### each row should only differ by vax data, including vaccine_weeks


# Exclusions --------------------------------------------------------------

# demog_labs_vax |>
#   filter(vax_series=="Primary series + 2 additional") |>
#   # filter(meets_bv_booster==T) |>
#   select(meets_bv_booster, vaccine_weeks) |>
#   table()

### meets_bv_booster (demog) aligns well with vaccine_weeks (vax)
### why some meets_bv_booster==F with vaccine_weeks <= 12? when filtered to "Primary series + 2 additional" ?
### these were not bivalent vaccines (bivalent==0 for all, must have been monovalent)
### answer below

# demog_labs_vax |>
#   filter(meets_bv_booster == F, vaccine_weeks <= 12) |>
#   filter(vax_series=="Primary series + 2 additional") |>
#   select(fake_mrn, vax_series, vaccine_weeks, bivalent) |>
#   summary()
### 2099 meets_bv_booster==F with vaccine_weeks <= 12? when filtered to "Primary series + 2 additional"

# demog_labs_vax |>
#   filter(vax_series=="Primary series + 2 additional") |>
#   filter(meets_bv_booster==T) |>
#   select(last_enroll_booster_weeks, vaccine_weeks) |>
#   mutate(diff=last_enroll_booster_weeks-vaccine_weeks) |>
#   summary()

### last_enroll_booster_weeks (demog) aligns well with vaccine_weeks (vax) when filtered to meets_bv_booster==T

fake_mrn_exclude_vax <- demog_labs_vax |>
  filter(meets_bv_booster == F, vaccine_weeks <= 12) |>
  filter(vax_series=="Primary series + 2 additional") |>
  select(fake_mrn)

# demog_labs_vax |> select(fake_mrn) |> unique() |> nrow() # 871285

# demog_labs_vax |>
#   filter(!fake_mrn %in% fake_mrn_exclude_vax$fake_mrn) |>
#   select(fake_mrn) |>
#   unique() |>
#   nrow() # 921779

# 923878-921779=2099

demog_labs_vax <- demog_labs_vax |>
  filter(!fake_mrn %in% fake_mrn_exclude_vax$fake_mrn)

# demog_labs_vax |> select(fake_mrn) |> unique() |> nrow() # 869275

# *** NEW check value moving forward 869275 ***


# Censoring Times ---------------------------------------------------------


### SAME for TWO and THREE level outcomes

### For A=0
### enr_end # end of KP enrollment IF < 78 weeks, or
### date of "Primary series + 2 additional" (tx devaition, 4 vax)

### For A=1
### enr_end # end of KP enrollment IF < 78 weeks, or
### date of "Primary series + 3 additional" (tx deviation, 5 vax)

### note for both A=0 and A=1

deviation_times <- demog_labs_vax |>
  filter(
    (meets_bv_booster == FALSE  & vax_series=="Primary series + 2 additional") |
      (meets_bv_booster == TRUE & vax_series=="Primary series + 3 additional")
  ) |>
  select(fake_mrn, vaccine_weeks) |>
  rename(deviation_t = vaccine_weeks)

# deviation_times$deviation_t |> summary()
### no NA values
# demog_labs_vax$enr_end_weeks |> summary()
### no NA values


##* retain for week matching
# Treatment times ---------------------------------------------------------


treatment_times <- demog_labs_vax |>
  filter(meets_bv_booster == TRUE & vax_series=="Primary series + 2 additional") |>
  select(fake_mrn, vaccine_weeks) |>
  rename(treatment_time = vaccine_weeks)


# Join Censoring Times ----------------------------------------------------


### note - cens time is actually per protocol censoring time,
###        itt censoring = enr_end_weeks

analysis_c <- demog_labs_vax |>
  left_join(deviation_times, by="fake_mrn") |>
  # take min of deviation_t and enr_end_weeks for each row to create censoring time
  # na.rm = T, else individuals without tx deviation times will be assigned NA
  mutate(cens_time = pmin(deviation_t, enr_end_weeks, na.rm = T)) |>
  # for each fake_mrn, all rows should be identical except for vaccination information
  # take 1 row per fake_mrn
  group_by(fake_mrn) |>
  slice_max(vax_num, with_ties = F) |>
  ungroup()

# analysis_c |> select(fake_mrn) |> unique() |> nrow() # 869275
# summary(analysis_c$cens_time)

### cross check with KP variable first_ne_dose_weeks

# analysis_c |>
#   mutate(diff=deviation_t-first_ne_dose_weeks) |>
#   filter(diff!=0) # |> nrow()

### 0 where deviation_t != first_ne_dose_weeks (after filtering above)
### *** exclusions done in section (Join Demog + Labs + Vax)


# Join Treatment Times ----------------------------------------------------


analysis_c_a <- analysis_c |>
  left_join(treatment_times, by="fake_mrn")

# analysis_c_a |> select(fake_mrn) |> unique() |> nrow() # 869275

# analysis_c_a |> filter(!is.na(treatment_time)) |> nrow()
# analysis_c_a |> filter(meets_bv_booster==T) |> nrow()
# matches


# Outcome Variables and Times ---------------------------------------------


### TWO level
# Y_TWO
# 0 = censor, 1 = test pos
# Y_TWO_t
# pmin(cens_time, lab_weeks)

table(analysis_c_a$result_2L)

### THREE level
# Y_THREE
# 0 = censor, 1 = test neg, 2 = test pos
# pmin(cens_time, lab_weeks)

table(analysis_c_a$result_3L)


### Censoring / end of FU distinction

### Y_TWO == 1 pos
### Y_TWO == 0 cens
### Y_TWO == 0 cens

### Y_THREE == 2 pos
### Y_THREE == 1 neg
### Y_THREE == 0 cens
### Y_THREE == 0 cens

### if enr_end_weeks == 78, then Y = 0
### if enr_end_weeks < 78, then Y = 0


### All times here are in CALENDAR time (weeks from 09/01)
##* Will be adjusted to index_time after matching


analysis_c_a_y <- analysis_c_a |>
  
  # two level outcome
  mutate(
    Y_TWO = case_when(
      result_2L == "positive" ~ 1,
      # is.na(result_2L) & pmin(cens_time, lab_weeks_2L, na.rm = TRUE) < 78 ~ 0,
      # is.na(result_2L) & pmin(cens_time, lab_weeks_2L, na.rm = TRUE) == 78 ~ NA_real_
      TRUE ~ 0
    ),
    Y_TWO_time_pp = pmin(cens_time, lab_weeks_2L, na.rm = TRUE), # cens at ltfu or tx dev
    Y_TWO_time_itt= pmin(enr_end_weeks, lab_weeks_2L, na.rm = TRUE) # cens at ltfu only
  ) |>
  
  # three level outcome
  mutate(
    Y_THREE = case_when(
      result_3L == "positive" ~ 2,
      result_3L == "negative" ~ 1,
      TRUE ~ 0
    ),
    Y_THREE_time_pp = pmin(cens_time, lab_weeks_3L, na.rm = TRUE), # cens at ltfu or tx dev
    Y_THREE_time_itt = pmin(enr_end_weeks, lab_weeks_3L, na.rm = TRUE) # cens at ltfu only
  )


# analysis_c_a_y |>
#   select(Y_TWO, Y_TWO_time_itt, Y_TWO_time_pp, Y_THREE, Y_THREE_time_itt, Y_THREE_time_pp) |>
#   summary()
# 
# analysis_c_a_y |> select(Y_TWO) |> table()
# analysis_c_a_y |> select(Y_THREE) |> table()


# Omit unnecessary variables ----------------------------------------------

# analysis_c_a_y |> summary()

analysis_data <- analysis_c_a_y |>
  select(-charlson_cat, # keep factor version
         -last_pe_booster, # same for all
         -enroll_doses, # unnecessary, 0 or 1
         -last_enroll_booster, # redundant with meets_bv_booster
         -has_bv_booster, # same as meets_bv_booster
         -infect, # using lab outcomes
         -months_last_vax_infect, # using last_vax_infect_weeks
         -flu_weeks, # colinear with flu_vax
         -last_pe_booster_weeks, # using last_vax_infect_weeks
         -last_enroll_booster_weeks, # don't adjust for this
         -first_ne_dose_weeks, # defined differently using deviation times above
         -covid_weeks, # using lab results
         -prior_inf_weeks, # just use last_vax_infect_weeks and adjust for last immun exposure
         -enr_start_weeks, # dont need to model
         -vax_name, # not adjust for make
         -vax_num, # don't need to model. slice_max earlier so just takes everyone's max value anyway
         -bivalent, # not informative after slice_max
         -vax_series, # not informative after slice_max
         -vaccine_weeks, # not informative after slice_max
         -deviation_t, # not informative after slice_max
         -lab_weeks_3L, # built into Y_X_time
         -result_3L, # built into Y_X
         -lab_weeks_2L, # built into Y_X_time
         -result_2L # built into Y_X
  ) |>
  mutate(treatment = case_when(meets_bv_booster == F ~ 0, meets_bv_booster == T ~ 1)) |>
  select(-meets_bv_booster)

# analysis_data |> summary()

##* treatment_time is retained: non-NA for A=1, NA for A=0
##* cens_time is retained: needed for risk-set eligibility

#netrisk#
write_rds(analysis_data, "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch.4/data_beforematch.rds")


# Weekly matching ---------------------------------------------------------

##* new process 2026-03-01

### At each vaccination week t during enrollment:
###   Treated  = individuals vaccinated during week t (floor(treatment_time) == t)
###   Controls = individuals NOT YET vaccinated at week t:
###              - never vaccinated (treatment_time is NA), OR
###              - vaccinated at a LATER week (floor(treatment_time) > t)
###            AND still enrolled in KP, event-free, and not already matched.
###
### Each individual appears EXACTLY ONCE in the final dataset.
### Once matched (as treated or control), they are removed from all future pools.
### A future-vaccinated individual matched as a control stays a control;
### they do not later enter as treated.
###
### Treatment is fixed at matching -- not time-varying over follow-up.
### For PP: future-vaccinated controls are censored at their vaccination time.
### For ITT: no censoring at crossover.
###
### NOTE: weeks are 0-indexed (week 0 = first week of enrollment starting 09/01)


### round treatment times to integer weeks for grouping
analysis_data <- analysis_data |>
  mutate(treat_week = if_else(treatment == 1, floor(treatment_time), NA_real_))

# verify distribution
# table(analysis_data$treat_week[analysis_data$treatment == 1])

vax_weeks <- sort(unique(na.omit(analysis_data$treat_week)))

cat("Beginning enrollment-period matching across weeks:", paste(vax_weeks, collapse=", "), "\n")

matched_pairs <- list()
already_matched <- character(0)  # tracks ALL matched individuals (treated + control)

for (t in vax_weeks) {
  
  # Treated: vaccinated during week t, not already matched as a control in an earlier week
  treated_t <- analysis_data |>
    filter(treat_week == t,
           enr_end_weeks > t,
           cens_time > t,                          #netrisk#
           Y_TWO_time_itt > t,                     #netrisk# was Y_THREE_time_itt in .3
           !fake_mrn %in% already_matched) |>
    mutate(treatment = 1L)
  
  # Eligible controls at week t:
  # 1. Not yet vaccinated: never vaccinated OR vaccinated at a later week
  # 2. Still enrolled AND no treatment deviation by week t
  # 3. No POSITIVE test before week t     #netrisk#
  # 4. Not already matched in an earlier week
  controls_t <- analysis_data |>
    filter(is.na(treat_week) | treat_week > t,   # not yet vaccinated at week t
           enr_end_weeks > t,                      # still enrolled in KP at week t
           cens_time > t,                          #netrisk#
           Y_TWO_time_itt > t,                     #netrisk# was Y_THREE_time_itt in .3 #netrisk#
           !fake_mrn %in% already_matched) |>
    mutate(treatment = 0L)
  
  cat("Week", t, ": ", nrow(treated_t), "treated,", nrow(controls_t), "eligible controls\n")
  
  if (nrow(treated_t) == 0 || nrow(controls_t) < 10) {
    cat("  Skipping week", t, "- insufficient pool\n")
    next
  }
  
  match_data_t <- bind_rows(treated_t, controls_t)
  
  m_t <- tryCatch(
    matchit(
      treatment ~ age_years + ndi + bmi + tests_count + last_vax_infect_weeks,
      data = match_data_t,
      method = "nearest",
      exact = ~ sex_admin + race + service_region + prior_inf + charlson_cat_fac # flu_vax removed 2026-07-28
    ),
    error = function(e) {
      cat("  matchit error at week", t, ":", conditionMessage(e), "\n")
      NULL
    }
  )
  
  if (is.null(m_t)) next
  
  md_t <- match.data(m_t) |>
    mutate(index_time = t)
  
  n_matched <- sum(md_t$treatment == 1)
  cat("  Matched", n_matched, "pairs at week", t, "\n")
  
  matched_pairs[[length(matched_pairs) + 1]] <- md_t
  
  # Mark ALL matched individuals as used -- removed from ALL future pools
  already_matched <- c(
    already_matched,
    md_t |> pull(fake_mrn)
  )
}

analysis_data_matched <- bind_rows(matched_pairs)

cat("\nTotal matched pairs:", sum(analysis_data_matched$treatment == 1), "\n") # Total matched pairs: 196881 

analysis_data_matched <- analysis_data_matched |>
  mutate(subclass = paste0("w", index_time, "_", subclass))

# verify all pairs are size 2
# analysis_data_matched |> count(subclass) |> count(n) # verified

# verify no individual appears more than once
# analysis_data_matched |> count(fake_mrn) |> filter(n > 1) |> nrow()  # = 0, verified


#netrisk#
# Neg test history for net risk -------------------------------------------

### neg_hist  = OPTION B input (recurrent negative intensity on {T2 > t})
### first_neg = OPTION A input (single-event NCO on its own risk set)

neg_hist <- analysis_data_matched |>
  select(fake_mrn, index_time, enr_end_weeks) |>
  inner_join(labs_neg_all, by = "fake_mrn") |>
  filter(neg_week >  index_time,        # after index only
         neg_week <= enr_end_weeks) |>  # inside enrollment; 
  mutate(neg_t = neg_week - index_time) |>
  select(fake_mrn, neg_t) |>
  arrange(fake_mrn, neg_t)

### Option A: time to FIRST negative after index.
### Positives are NON-truncating here -- someone who tests positive at week 6
### and never tests negative contributes panel-N rows all the way to censoring.
first_neg <- neg_hist |>
  group_by(fake_mrn) |>
  slice_min(neg_t, with_ties = FALSE) |>
  ungroup() |>
  rename(neg1_t = neg_t)

### CHECK: multi-negative weeks decide binary vs Poisson for Option B.
### If mostly n == 1, binary N1 loses nothing and keeps both stages
### in the same model family (binomial). If not, use Poisson + log-exposure
### offset for stage one, and note the EQC log-rate/log-odds mismatch.
# neg_hist |> count(fake_mrn, wk = ceiling(neg_t)) |> count(n > 1) mostly n=1


# Adjust times to index time after matching -------------------------------

analysis_data_matched_adj <- analysis_data_matched |>
  mutate(
    
    # --- ITT: no censoring at vaccination crossover ---
    Y_TWO_time_itt_adj   = Y_TWO_time_itt - index_time,
    Y_THREE_time_itt_adj = Y_THREE_time_itt - index_time,
    # outcomes unchanged for ITT
    Y_TWO_itt   = Y_TWO,
    Y_THREE_itt = Y_THREE,
    
    # --- PP: censor future-vaccinated controls at their vaccination time ---
    # For controls with non-NA treatment_time, vaccination = treatment deviation
    # If treatment_time < original PP event/cens time, cap time and censor event
    # handles censoring during enrollment period
    # censoring after enrollment period handled upstream
    
    control_vax_cens = case_when(
      treatment == 0 & !is.na(treatment_time) ~ treatment_time,
      TRUE ~ Inf
    ),
    
    C_itt = enr_end_weeks - index_time, #netrisk# censoring horizon
    C_pp  = pmin(cens_time, control_vax_cens) - index_time, #netrisk# censoring horizon
    
    # Cap PP times at vaccination for future-vaccinated controls
    Y_TWO_time_pp_adj   = pmin(Y_TWO_time_pp, control_vax_cens) - index_time,
    Y_THREE_time_pp_adj = pmin(Y_THREE_time_pp, control_vax_cens) - index_time,
    
    # If PP time was bound by deviation/admin censoring before the event,
    # censor the event. Covers BOTH arms (treated P+3, control P+2),
    # inside or outside enrollment window.
    Y_TWO_pp = case_when(
      Y_TWO == 1L & Y_TWO_time_pp < Y_TWO_time_itt ~ 0L,
      treatment == 0 & !is.na(treatment_time) & treatment_time < Y_TWO_time_pp ~ 0L,
      TRUE ~ Y_TWO
    ),
    Y_THREE_pp = case_when(
      Y_THREE != 0L & Y_THREE_time_pp < Y_THREE_time_itt ~ 0L,
      treatment == 0 & !is.na(treatment_time) & treatment_time < Y_THREE_time_pp ~ 0L,
      TRUE ~ Y_THREE
    )
    
  ) |>
  select(-control_vax_cens)

# Sanity checks:
# analysis_data_matched_adj |> filter(Y_TWO_time_itt_adj <= 0) |> nrow()   # 0
# analysis_data_matched_adj |> filter(Y_THREE_time_itt_adj <= 0) |> nrow() # 14594
# analysis_data_matched_adj |> filter(Y_TWO_time_pp_adj < 0) |> nrow()     # 0
# analysis_data_matched_adj |> filter(Y_THREE_time_pp_adj < 0) |> nrow() # 12624
# 
# # How many future-vaccinated controls had their events censored?
# analysis_data_matched_adj |> filter(treatment == 0, Y_TWO != Y_TWO_pp) |> nrow() # 3385 #netrisk#
# analysis_data_matched_adj |> filter(treatment == 0, Y_THREE != Y_THREE_pp) |> nrow() # 9616 #netrisk#

##* Handle remaining negative/zero adjusted times 
##* This should be rare given risk-set eligibility filters.
##* Drop affected pairs to be safe.

analysis_data_matched_adj <- analysis_data_matched_adj |>
  group_by(subclass) |>
  mutate(pair_valid = all(Y_TWO_time_itt_adj > 0) &
           all(Y_TWO_time_pp_adj  > 0)) |>   #netrisk# Y_THREE conditions removed
  ungroup()

cat("Pairs dropped due to non-positive adjusted times:",
    sum(!analysis_data_matched_adj$pair_valid) / 2, "\n") # 0

analysis_data_matched_adj <- analysis_data_matched_adj |>
  filter(pair_valid) |>
  select(-pair_valid)

#netrisk#
# Drop negatives from dropped pairs -----------------------------------------

# Keep neg_hist consistent with the final matched cohort.
neg_hist <- neg_hist |>
  filter(fake_mrn %in% analysis_data_matched_adj$fake_mrn)

# Split into Y2 and Y3 datasets -------------------------------------------


data_Y2 <- analysis_data_matched_adj |>
  select(-Y_THREE, -Y_THREE_itt, -Y_THREE_pp,
         -Y_THREE_time_itt, -Y_THREE_time_pp,
         -Y_THREE_time_itt_adj, -Y_THREE_time_pp_adj) |>
  rename(Y2     = Y_TWO_itt,
         Y2_pp  = Y_TWO_pp,
         Y2_itt_t = Y_TWO_time_itt_adj,
         Y2_pp_t = Y_TWO_time_pp_adj,
         Y2_itt_t_unadj = Y_TWO_time_itt,
         Y2_pp_t_unadj = Y_TWO_time_pp) |>
  select(-Y_TWO)

  #netrisk# data_Y2 is UNCHANGED in structure -- it is already the net-risk positive
  #netrisk# panel span, min(T2, C), and is reused by all four analyses
  #netrisk# (EQC/PCI x Option A/B) as the stage-two dataset.

  #netrisk# data_Y3 NOT built here. The Y_TWO_time_itt eligibility change admits
  #netrisk# individuals with a pre-index negative test, for whom
  #netrisk# Y_THREE_time_itt_adj <= 0 (14,594 individuals, no longer removed by
  #netrisk# pair_valid). They are valid under net risk but were correctly
  #netrisk# ineligible under competing risks.
  #netrisk# CR comparison arm uses data_weekmatch.3/data_Y3_weekmatch.rds

#netrisk#
# data_OA: OPTION A panel N -----------------------------------------------

### Single-event NCO: time to first negative after index, on its OWN risk
### set {T1^A > t}. A positive test is NOT an event and does NOT truncate.

data_OA <- analysis_data_matched_adj |>
  select(-Y_TWO, -Y_TWO_itt, -Y_TWO_pp,
         -Y_TWO_time_itt, -Y_TWO_time_pp,
         -Y_TWO_time_itt_adj, -Y_TWO_time_pp_adj,
         -Y_THREE, -Y_THREE_itt, -Y_THREE_pp,
         -Y_THREE_time_itt, -Y_THREE_time_pp,
         -Y_THREE_time_itt_adj, -Y_THREE_time_pp_adj) |>
  left_join(first_neg, by = "fake_mrn") |>
  mutate(
    # ITT: censor at end of enrollment only
    Y1_itt_t = pmin(neg1_t, C_itt, na.rm = TRUE),
    Y1_itt   = if_else(!is.na(neg1_t) & neg1_t <= C_itt, 1L, 0L),
    # PP: censor at deviation / crossover as well
    Y1_pp_t  = pmin(neg1_t, C_pp, na.rm = TRUE),
    Y1_pp    = if_else(!is.na(neg1_t) & neg1_t <= C_pp, 1L, 0L)
  )

# data_OA |> filter(Y1_itt_t <= 0) |> nrow()  # expect 0
# data_OA |> count(Y1_itt)
# data_OA |> count(Y1_pp)


# Truncate data for shorter follow-up -------------------------------------


max_follow <- 52

data_Y2 <- data_Y2 |>
  mutate(
    # itt truncation
    Y2_itt_t_trunc = pmin(Y2_itt_t, max_follow),
    Y2_itt_trunc   = if_else(Y2_itt_t > max_follow & Y2 == 1L, 0L, Y2),
    
    # pp truncation
    Y2_pp_t_trunc  = pmin(Y2_pp_t, max_follow),
    Y2_pp_trunc    = if_else(Y2_pp_t > max_follow & Y2_pp == 1L, 0L, Y2_pp)
  )


#netrisk#
data_OA <- data_OA |>
  mutate(
    # itt truncation
    Y1_itt_t_trunc = pmin(Y1_itt_t, max_follow),
    Y1_itt_trunc   = if_else(Y1_itt_t > max_follow & Y1_itt == 1L, 0L, Y1_itt),
    
    # pp truncation
    Y1_pp_t_trunc  = pmin(Y1_pp_t, max_follow),
    Y1_pp_trunc    = if_else(Y1_pp_t > max_follow & Y1_pp == 1L, 0L, Y1_pp)
  )

#netrisk#
neg_hist <- neg_hist |>
  filter(neg_t <= max_follow)


# Write Finalized Data ----------------------------------------------------

### first converting to data.table

setDT(data_Y2)
setDT(data_OA)     #netrisk#
setDT(neg_hist)    #netrisk#

write_csv(data_Y2, "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch.4/data_Y2_weekmatch.csv")
write_rds(data_Y2, "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch.4/data_Y2_weekmatch.rds")
write_csv(data_OA,  "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch.4/data_OA_weekmatch.csv")
write_rds(data_OA,  "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch.4/data_OA_weekmatch.rds")
write_csv(neg_hist, "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch.4/neg_hist_weekmatch.csv")
write_rds(neg_hist, "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch.4/neg_hist_weekmatch.rds")

### last written 2026-08-29 at 0903 (net risk additions)


# Clean Up Environment ----------------------------------------------------


rm(list = setdiff(ls(), c("data_Y2", "data_OA", "neg_hist")))  #netrisk#






