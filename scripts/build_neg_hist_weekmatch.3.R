##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Build neg_hist for the data_weekmatch.3 cohort
##----- Adds the Option B net-risk input to the existing CR cohort.
##----- No rematching: pairs, index weeks and eligibility are unchanged.
##----- last updated 2026-09-18


# Packages ----------------------------------------------------------------


library(tidyverse)
library(data.table)


# Data --------------------------------------------------------------------


labs <- read_csv("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/labs_250701.csv")

data_path_3 <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch.3/"

data_Y2_3 <- read_rds(paste0(data_path_3, "data_Y2_weekmatch.rds"))

### data_Y2 is the final matched cohort for .3, so joining to it applies the
### pair_valid drops without needing to rebuild analysis_data_matched_adj.
stopifnot(all(c("fake_mrn", "index_time", "enr_end_weeks") %in% names(data_Y2_3)))
stopifnot(!anyDuplicated(data_Y2_3$fake_mrn))


# Labs / Testing ----------------------------------------------------------


### identical to data_cleaning_weekmatch_net.R

labs_clean <- labs |>
  
  filter(test_type=="SARS_COV_2_NAAT",
         pt_loc=="O",
         result!="no result") |>
  
  filter(lab_weeks >= 0) # tests eligible as outcomes need positive value

# retain all negative tests for net risk
labs_neg_all <- labs_clean |>
  filter(result == "negative") |>
  select(fake_mrn, lab_weeks) |>
  rename(neg_week = lab_weeks)


# neg_hist ----------------------------------------------------------------


### neg_hist = OPTION B input (recurrent negative intensity on {T2 > t})

neg_hist <- data_Y2_3 |>
  select(fake_mrn, index_time, enr_end_weeks) |>
  inner_join(labs_neg_all, by = "fake_mrn") |>
  filter(neg_week >  index_time,        # after index only
         neg_week <= enr_end_weeks) |>  # inside enrollment
  mutate(neg_t = neg_week - index_time) |>
  select(fake_mrn, neg_t) |>
  arrange(fake_mrn, neg_t)


# Truncate data for shorter follow-up -------------------------------------


max_follow <- 52

neg_hist <- neg_hist |>
  filter(neg_t <= max_follow)


# Write -------------------------------------------------------------------


setDT(neg_hist)

write_csv(neg_hist, paste0(data_path_3, "neg_hist_weekmatch.csv"))
write_rds(neg_hist, paste0(data_path_3, "neg_hist_weekmatch.rds"))
