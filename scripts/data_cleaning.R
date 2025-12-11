##----- Beau Schaeffer
##----- Kaiser TTE-TND
##----- Data cleaning

# Packages ----------------------------------------------------------------

library(tidyverse)
# install.packages("MatchIt")
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
    # significant number, imputation later *****
    
    last_enroll_booster=="NONE" | last_enroll_booster=="Primary series + 2 additional",
    # keep individuals who were P+1 or P+2 at the time of their last enrollment booster 
    # I believe "NONE" indicates P+1
    # 949070 (894937 after charlson filter)
    
    months_last_vax_infect >= 3.0,
    #  ????? last immunological event at least 3 months prior to enroll
    # redundant with last_pe_booster_weeks?
    # 923878 (871285 after charlson filter)
    
    last_pe_booster_weeks < -12.0
    # ????? last booster at least 3 months prior to enroll
    # redundant with months_last_vax_infect?
    # 923878 (871285 after charlson filter)
    
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

demog_labs <- demog_clean |> 
  left_join(labs_clean_3L, by = "fake_mrn") |> 
  left_join(labs_clean_2L, by = "fake_mrn")

  # demog_labs |> select(fake_mrn) |> unique() |> nrow() # 923878 (871285 after charlston filter)


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

  # demog_labs_vax |> select(fake_mrn) |> unique() |> nrow() # 923878 (871285 after charlston filter)

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

  # demog_labs_vax |> select(fake_mrn) |> unique() |> nrow() # 923878 (871285 after charlson)

  # demog_labs_vax |> 
  #   filter(!fake_mrn %in% fake_mrn_exclude_vax$fake_mrn) |> 
  #   select(fake_mrn) |> 
  #   unique() |> 
  #   nrow() # 921779
  
  # 923878-921779=2099

demog_labs_vax <- demog_labs_vax |> 
  filter(!fake_mrn %in% fake_mrn_exclude_vax$fake_mrn)

  # demog_labs_vax |> select(fake_mrn) |> unique() |> nrow() # 921779 (869275 after charlson)

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


# Event Time Corrections --------------------------------------------------


  # for A=1, need to retain time of vaccination
  # will need to subtract this difference after matching so that both pairs times start from A=1 time zero

treatment_times <- demog_labs_vax |>
  filter(meets_bv_booster == TRUE & vax_series=="Primary series + 2 additional") |>
  select(fake_mrn, vaccine_weeks) |> 
  rename(treatment_time = vaccine_weeks)

  # considerations
  # treatment times range from 0-12
  # only want to match to unvaccinated who are still eligible at time of vaccination
  # A=0 censoring times range from 4-78 weeks
  # analysis_c_y |> filter(meets_bv_booster==FALSE) |> select(cens_time) |> filter(cens_time<=12) |> nrow()
  # 14450/921779 = 0.01567621 --> proportion A=0 censored before 12 weeks
  # let's just match and exclude pairs where vaccination time for A=1 > censoring time A=0


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


# Join Event Time Corrections ---------------------------------------------

analysis_c_a <- analysis_c |> 
  left_join(treatment_times, by="fake_mrn")

  # analysis_c_a |> select(fake_mrn) |> unique() |> nrow() # 921779 (869275 after charlson)


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


  analysis_c_a_y |>
    select(Y_TWO, Y_TWO_time_itt, Y_TWO_time_pp, Y_THREE, Y_THREE_time_itt, Y_THREE_time_pp) |>
    summary()
  
  analysis_c_a_y |> select(Y_TWO) |> table()
  analysis_c_a_y |> select(Y_THREE) |> table()


# Omit unnecessary variables ----------------------------------------------

analysis_c_a_y |> summary()

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
         -enr_end_weeks, # dont need to model
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

analysis_data |> summary()


# Matching ----------------------------------------------------------------


m.out <- matchit(
  treatment ~ age_years + ndi + bmi + tests_count + last_vax_infect_weeks,
  data = analysis_data,
  method = "nearest",
  exact = ~ sex_admin + race + service_region + flu_vax + prior_inf + charlson_cat_fac
)

summary(m.out)

analysis_data_matched <- match.data(m.out)


  # make sure all groups size 2
  # analysis_data_matched |> select(subclass, fake_mrn) |> group_by(subclass) |> mutate(n=n()) |> summary()


# Adjusting times after matching ------------------------------------------


  # Y_TWO_time and Y_THREE_time
  # for A=0, time from first eligible day to event
  # for A=1, time from vaccination to event
  # for pairs, will need to adjust time for A=0, A=1 should be correct

  # analysis_data_matched |>
  #   select(fake_mrn, subclass, treatment, treatment_time,
  #          Y_TWO, Y_TWO_time_itt, Y_TWO_time_pp, Y_THREE, Y_THREE_time_itt, Y_THREE_time_pp) |>
  #   arrange(subclass) |>
  #   print(n=100)


analysis_data_matched_adj <- analysis_data_matched |> 
  group_by(subclass) |> 
  
  # select(fake_mrn, subclass, treatment, treatment_time, Y_TWO, Y_TWO_time, Y_THREE, Y_THREE_time) |> # TOGGLE
  
  # 2025-12-10 - realize should be making this adjustment to both A=0 and A=1, not just A=0
  # mutate(subclass_treat_time = max(treatment_time[treatment == 1], na.rm = TRUE),
  #        Y_TWO_time_itt_adj = if_else(treatment == 0, Y_TWO_time_itt - subclass_treat_time, Y_TWO_time_itt),
  #        Y_TWO_time_pp_adj = if_else(treatment == 0, Y_TWO_time_pp - subclass_treat_time, Y_TWO_time_pp),
  #        Y_THREE_time_itt_adj = if_else(treatment == 0, Y_THREE_time_itt - subclass_treat_time, Y_THREE_time_itt),
  #        Y_THREE_time_pp_adj = if_else(treatment == 0, Y_THREE_time_pp - subclass_treat_time, Y_THREE_time_pp)) |> 
  
  # 2025-12-10 - editing above section to adjust follow-up time for both members of matched pairs
  mutate(subclass_treat_time = max(treatment_time[treatment == 1], na.rm = TRUE),
         Y_TWO_time_itt_adj = Y_TWO_time_itt - subclass_treat_time,
         Y_TWO_time_pp_adj = Y_TWO_time_pp - subclass_treat_time,
         Y_THREE_time_itt_adj = Y_THREE_time_itt - subclass_treat_time,
         Y_THREE_time_pp_adj = Y_THREE_time_pp - subclass_treat_time) |> 
  
  ungroup() |> 
  select(-subclass_treat_time)


  analysis_data_matched_adj |>
    filter(Y_TWO_time_pp_adj < 0 | Y_THREE_time_pp_adj < 0) |>
    # filter(Y_TWO_time_itt_adj <= 0 | Y_THREE_time_itt_adj <= 0) |>
    select(subclass) |>
    unique() |>
    nrow() # 11762 (11690 after charlston) (20806 after additional time adjustment)

  # analysis_data_matched_adj |> select(subclass) |> unique() |> nrow() # 246679 (236276)

  # 11762/246679=0.0476814
  # 11690/236276=0.04947604
  # 20806/236276=0.08805803

  # 8.8% of matches where A=1 matched to ineligible A=0 control
  # censor these pairs

analysis_data_Y_TWO <- analysis_data_matched_adj |> 
  select(-Y_THREE, -Y_THREE_time_itt, -Y_THREE_time_pp, -Y_THREE_time_itt_adj, -Y_THREE_time_pp_adj) |> 
  mutate(valid_pair = if_else(Y_TWO_time_itt_adj < 0, FALSE, TRUE)) |>
  group_by(subclass) |> 
  mutate(pair_invalid = any(!valid_pair)) |> 
  ungroup() |> 
  filter(!pair_invalid) |>
  select(-valid_pair, -pair_invalid)

table(analysis_data_matched_adj$subclass %in% analysis_data_Y_TWO$subclass) # 11120 censored

analysis_data_Y_THREE <- analysis_data_matched_adj |> 
  select(-Y_TWO, -Y_TWO_time_itt, -Y_TWO_time_pp, -Y_TWO_time_itt_adj, -Y_TWO_time_pp_adj) |> 
  mutate(valid_pair = if_else(Y_THREE_time_itt_adj < 0, FALSE, TRUE)) |> 
  group_by(subclass) |> 
  mutate(pair_invalid = any(!valid_pair)) |> 
  ungroup() |> 
  filter(!pair_invalid) |> 
  select(-valid_pair, -pair_invalid)

table(analysis_data_matched_adj$subclass %in% analysis_data_Y_THREE$subclass) # 41612 censored


# Clean Up Variable Names -------------------------------------------------

names(analysis_data_Y_TWO)

data_Y2 <- analysis_data_Y_TWO |> 
  rename(Y2 = Y_TWO,
         Y2_itt_t = Y_TWO_time_itt_adj,
         Y2_pp_t = Y_TWO_time_pp_adj,
         Y2_itt_t_unadj = Y_TWO_time_itt,
         Y2_pp_t_unadj = Y_TWO_time_pp)


names(analysis_data_Y_THREE)

data_Y3 <- analysis_data_Y_THREE |> 
  rename(Y3 = Y_THREE,
         Y3_itt_t = Y_THREE_time_itt_adj,
         Y3_pp_t = Y_THREE_time_pp_adj,
         Y3_itt_t_unadj = Y_THREE_time_itt,
         Y3_pp_t_unadj = Y_THREE_time_pp)


# Truncate data for shorter follow-up -------------------------------------

max_follow <- 52

data_Y2 <- data_Y2 |>
  mutate(
    # itt truncation
    Y2_itt_t_trunc = pmin(Y2_itt_t, max_follow),
    Y2_itt_trunc   = if_else(Y2_itt_t > max_follow & Y2 == 1L, 0L, Y2),
    
    # pp truncation
    Y2_pp_t_trunc  = pmin(Y2_pp_t, max_follow),
    Y2_pp_trunc    = if_else(Y2_pp_t > max_follow & Y2 == 1L, 0L, Y2)
  )


data_Y3 <- data_Y3 |>
  mutate(
    # itt truncation
    Y3_itt_t_trunc = pmin(Y3_itt_t, max_follow),
    Y3_itt_trunc   = if_else(Y3_itt_t > max_follow & Y3 != 0L, 0L, Y3),
    
    # pp truncation
    Y3_pp_t_trunc  = pmin(Y3_pp_t, max_follow),
    Y3_pp_trunc    = if_else(Y3_pp_t > max_follow & Y3 != 0L, 0L, Y3)
  )



# Write Finalized Data ----------------------------------------------------

  ### first converting to data.table

# class(analysis_data_Y_TWO)
# format(object.size(analysis_data_Y_TWO), units = "auto")
setDT(data_Y3)
# class(analysis_data_Y_TWO)
# format(object.size(analysis_data_Y_TWO), units = "auto")

# class(analysis_data_Y_THREE)
# format(object.size(analysis_data_Y_THREE), units = "auto")
setDT(data_Y2)
# class(analysis_data_Y_THREE)
# format(object.size(analysis_data_Y_THREE), units = "auto")

# write_csv(data_Y2, "cleaned_data/data_Y2.csv")
# write_rds(data_Y2, "cleaned_data/data_Y2.rds")
# 
# write_csv(data_Y3, "cleaned_data/data_Y3.csv")
# write_rds(data_Y3, "cleaned_data/data_Y3.rds")

  ### last written 2025-07-25 at 1425
  ### last written 2025-07-29 at 1630 (charlston filter, data.frame)
  ### last written 2025-07-30 at 1608 (corrected censoring indicator)
  ### last written 2025-07-30 at 1638 (indicator didnt need changing after all)
  ### last written 2025-08-14 at 1153 (distinguishing itt and pp effects, event and censoring times)
  ### last written 2025-12-05 at 1543 (shorten FU time to 52 weeks)
  ### last written 2025-12-05 at 1610 (revert back)
  ### last written 2025-12-08 at 1425 (shorten FU time to 52 weeks, data good, was a model issue)
  ### last written 2025-12-09 at 1457 (properly shorten FU time to 52 weeks)
  ### last written 2025-12-10 at 1551 (adjust follow-up time of A=1 in pairs; overlooked earlier)



# Clean Up Environment ----------------------------------------------------


rm(list = setdiff(ls(), c("data_Y2", "data_Y3")))






