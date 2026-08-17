##----- Beau Schaeffer
##----- Kaiser TTE-TND
##----- Covariate balance tables

library(tidyverse)
library(tableone)
library(kableExtra)


# Load data ---------------------------------------------------------------


data_beforematch <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch.3/data_beforematch.rds")
data_Y3 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch.3/data_Y3_weekmatch.rds")


# Define vars -------------------------------------------------------------


vars <- c("age_years", "bmi", "ndi", "tests_count", "last_vax_infect_weeks",
          "sex_admin", "race", "service_region", "flu_vax",
          "prior_inf", "charlson_cat_fac")

factor_vars <- c("sex_admin", "race", "service_region", "flu_vax",
                 "prior_inf", "charlson_cat_fac")

# weeks <- sort(unique(data_Y3$index_time))


# Covariate balance before matching ---------------------------------------


analysis_data <- data_beforematch |>
  mutate(treat_week = if_else(treatment == 1, floor(treatment_time), NA_real_))

matched_ledger <- data_Y3 |>
  select(fake_mrn, index_time) |>
  distinct()

vax_weeks <- sort(unique(na.omit(analysis_data$treat_week)))

ever_tx_pool <- character(0)   # ever entered a treated pool
ever_ct_pool <- character(0)   # ever entered a control pool
pool_log     <- list()

for (t in vax_weeks) {
  
  already <- matched_ledger$fake_mrn[matched_ledger$index_time < t]
  
  elig_t <- analysis_data |>
    filter(enr_end_weeks   > t,
           Y_THREE_time_itt > t,
           !fake_mrn %in% already)
  
  tx_t <- elig_t |> filter(treat_week == t) |> pull(fake_mrn)
  ct_t <- elig_t |> filter(is.na(treat_week) | treat_week > t) |> pull(fake_mrn)
  
  ever_tx_pool <- union(ever_tx_pool, tx_t)
  ever_ct_pool <- union(ever_ct_pool, ct_t)
  
  pool_log[[as.character(t)]] <- tibble(week = t,
                                        n_treated  = length(tx_t),
                                        n_controls = length(ct_t))
  
  cat("Week", t, ": ", length(tx_t), "treated,", length(ct_t), "eligible controls\n")
}

pool_log <- bind_rows(pool_log)

data_beforematch_pool <- analysis_data |>
  filter(fake_mrn %in% c(ever_tx_pool, ever_ct_pool)) |>
  mutate(
    pool_role = if_else(fake_mrn %in% ever_tx_pool, 1L, 0L),
    treatment_text = factor(pool_role, levels = c(1, 0),
                            labels = c("Booster", "No Booster"))
  )

# fate decomposition

data_beforematch |>
  filter(treatment == 1) |>
  mutate(treat_week = floor(treatment_time),
         fate = case_when(
           fake_mrn %in% matched_ledger$fake_mrn[
             matched_ledger$fake_mrn %in% ever_tx_pool]      ~ "matched as treated",
           !fake_mrn %in% ever_tx_pool &
             fake_mrn %in% ever_ct_pool                      ~ "consumed as control",
           Y_THREE_time_itt <= treat_week                    ~ "tested before vax week",
           enr_end_weeks    <= treat_week                    ~ "disenrolled before vax week",
           TRUE                                             ~ "other"
         )) |>
  count(fate)


# Covariate balance before matching ---------------------------------------


balance_before <- CreateTableOne(
  vars        = vars,
  factorVars  = factor_vars,
  strata      = "treatment_text",
  data        = data_beforematch_pool,
  addOverall  = FALSE
)

balance_mat_before <- print(balance_before, smd = TRUE, showAllLevels = TRUE, printToggle = FALSE, quote = FALSE, noSpaces = TRUE)
balance_mat_before <- balance_mat_before[, colnames(balance_mat_before) != "test"]

rownames(balance_mat_before) <- gsub("age_years (mean (SD))", "Age (mean (SD))", rownames(balance_mat_before), fixed = TRUE)
rownames(balance_mat_before) <- gsub("bmi (mean (SD))", "BMI (mean (SD))", rownames(balance_mat_before), fixed = TRUE)
rownames(balance_mat_before) <- gsub("ndi (mean (SD))", "NDI (mean (SD))", rownames(balance_mat_before), fixed = TRUE)
rownames(balance_mat_before) <- gsub("tests_count (mean (SD))", "Test count (mean (SD))", rownames(balance_mat_before), fixed = TRUE)
rownames(balance_mat_before) <- gsub("last_vax_infect_weeks (mean (SD))", "Weeks last inf/vax (mean (SD))", rownames(balance_mat_before), fixed = TRUE)
rownames(balance_mat_before) <- gsub("sex_admin (%)", "Sex (%)", rownames(balance_mat_before), fixed = TRUE)
rownames(balance_mat_before) <- gsub("race (%)", "Race (%)", rownames(balance_mat_before), fixed = TRUE)
rownames(balance_mat_before) <- gsub("service_region (%)", "Service region (%)", rownames(balance_mat_before), fixed = TRUE)
rownames(balance_mat_before) <- gsub("flu_vax (%)", "Prior year flu (%)", rownames(balance_mat_before), fixed = TRUE)
rownames(balance_mat_before) <- gsub("prior_inf (%)", "Prior COVID-19 infection (%)", rownames(balance_mat_before), fixed = TRUE)
rownames(balance_mat_before) <- gsub("charlson_cat_fac (%)", "Charlson (%)", rownames(balance_mat_before), fixed = TRUE)

var_names_before <- rownames(balance_mat_before)
linesep_vec_before <- rep("", nrow(balance_mat_before))
is_new_var_before <- c(TRUE, !grepl("^ ", var_names_before[-1]) & var_names_before[-1] != "")
is_last_of_group_before <- c(is_new_var_before[-1], TRUE)
linesep_vec_before[is_last_of_group_before] <- "\\addlinespace"

kbl_out_before <- kable(balance_mat_before, format = "latex", booktabs = TRUE, escape = TRUE, linesep = linesep_vec_before,
                         caption = "Covariate balance before matching") |>
  kable_styling(latex_options = c("HOLD_position", "scale_down"))

kbl_str_before <- as.character(kbl_out_before)
kbl_str_before <- gsub("<0.001", "$<$0.001", kbl_str_before, fixed = TRUE)

cat(kbl_str_before)


# Covariate balance after matching ----------------------------------------


data_Y3 <- data_Y3 |>
  mutate(treatment_text = factor(treatment, levels = c(1, 0),
                                 labels = c("Booster", "No Booster")))


tbl_after <- CreateTableOne(
  vars        = vars,
  factorVars  = factor_vars,
  strata      = "treatment_text",
  data        = data_Y3,
  addOverall  = FALSE
)

tbl_mat_after <- print(tbl_after, smd = TRUE, showAllLevels = TRUE, printToggle = FALSE, quote = FALSE, noSpaces = TRUE)
tbl_mat_after <- tbl_mat_after[, colnames(tbl_mat_after) != "test"]

rownames(tbl_mat_after) <- gsub("age_years (mean (SD))", "Age (mean (SD))", rownames(tbl_mat_after), fixed = TRUE)
rownames(tbl_mat_after) <- gsub("bmi (mean (SD))", "BMI (mean (SD))", rownames(tbl_mat_after), fixed = TRUE)
rownames(tbl_mat_after) <- gsub("ndi (mean (SD))", "NDI (mean (SD))", rownames(tbl_mat_after), fixed = TRUE)
rownames(tbl_mat_after) <- gsub("tests_count (mean (SD))", "Test count (mean (SD))", rownames(tbl_mat_after), fixed = TRUE)
rownames(tbl_mat_after) <- gsub("last_vax_infect_weeks (mean (SD))", "Weeks last inf/vax (mean (SD))", rownames(tbl_mat_after), fixed = TRUE)
rownames(tbl_mat_after) <- gsub("sex_admin (%)", "Sex (%)", rownames(tbl_mat_after), fixed = TRUE)
rownames(tbl_mat_after) <- gsub("race (%)", "Race (%)", rownames(tbl_mat_after), fixed = TRUE)
rownames(tbl_mat_after) <- gsub("service_region (%)", "Service region (%)", rownames(tbl_mat_after), fixed = TRUE)
rownames(tbl_mat_after) <- gsub("flu_vax (%)", "Prior year flu (%)", rownames(tbl_mat_after), fixed = TRUE)
rownames(tbl_mat_after) <- gsub("prior_inf (%)", "Prior COVID-19 infection (%)", rownames(tbl_mat_after), fixed = TRUE)
rownames(tbl_mat_after) <- gsub("charlson_cat_fac (%)", "Charlson (%)", rownames(tbl_mat_after), fixed = TRUE)

var_names_after <- rownames(tbl_mat_after)
linesep_vec_after <- rep("", nrow(tbl_mat_after))
is_new_var_after <- c(TRUE, !grepl("^ ", var_names_after[-1]) & var_names_after[-1] != "")
is_last_of_group_after <- c(is_new_var_after[-1], TRUE)
linesep_vec_after[is_last_of_group_after] <- "\\addlinespace"

kbl_out_after <- kable(tbl_mat_after, format = "latex", booktabs = TRUE, escape = TRUE, linesep = linesep_vec_after,
                         caption = "Covariate balance after matching") |>
  kable_styling(latex_options = c("HOLD_position", "scale_down"))

kbl_str_after <- as.character(kbl_out_after)
kbl_str_after <- gsub("<0.001", "$<$0.001", kbl_str_after, fixed = TRUE)

cat(kbl_str_after)


# Covariate balance end of follow up --------------------------------------

max_follow <- 52

data_endfu_itt <- data_Y3 |> filter(Y3_itt_t > max_follow)
data_endfu_pp  <- data_Y3 |> filter(Y3_pp_t  > max_follow)

# itt

tbl_end_itt <- CreateTableOne(
  vars        = vars,
  factorVars  = factor_vars,
  strata      = "treatment_text",
  data        = data_endfu_itt,
  addOverall  = FALSE
)

tbl_mat_end_itt <- print(tbl_end_itt, smd = TRUE, showAllLevels = TRUE, printToggle = FALSE, quote = FALSE, noSpaces = TRUE)
tbl_mat_end_itt <- tbl_mat_end_itt[, colnames(tbl_mat_end_itt) != "test"]

rownames(tbl_mat_end_itt) <- gsub("age_years (mean (SD))", "Age (mean (SD))", rownames(tbl_mat_end_itt), fixed = TRUE)
rownames(tbl_mat_end_itt) <- gsub("bmi (mean (SD))", "BMI (mean (SD))", rownames(tbl_mat_end_itt), fixed = TRUE)
rownames(tbl_mat_end_itt) <- gsub("ndi (mean (SD))", "NDI (mean (SD))", rownames(tbl_mat_end_itt), fixed = TRUE)
rownames(tbl_mat_end_itt) <- gsub("tests_count (mean (SD))", "Test count (mean (SD))", rownames(tbl_mat_end_itt), fixed = TRUE)
rownames(tbl_mat_end_itt) <- gsub("last_vax_infect_weeks (mean (SD))", "Weeks last inf/vax (mean (SD))", rownames(tbl_mat_end_itt), fixed = TRUE)
rownames(tbl_mat_end_itt) <- gsub("sex_admin (%)", "Sex (%)", rownames(tbl_mat_end_itt), fixed = TRUE)
rownames(tbl_mat_end_itt) <- gsub("race (%)", "Race (%)", rownames(tbl_mat_end_itt), fixed = TRUE)
rownames(tbl_mat_end_itt) <- gsub("service_region (%)", "Service region (%)", rownames(tbl_mat_end_itt), fixed = TRUE)
rownames(tbl_mat_end_itt) <- gsub("flu_vax (%)", "Prior year flu (%)", rownames(tbl_mat_end_itt), fixed = TRUE)
rownames(tbl_mat_end_itt) <- gsub("prior_inf (%)", "Prior COVID-19 infection (%)", rownames(tbl_mat_end_itt), fixed = TRUE)
rownames(tbl_mat_end_itt) <- gsub("charlson_cat_fac (%)", "Charlson (%)", rownames(tbl_mat_end_itt), fixed = TRUE)

var_names_end_itt <- rownames(tbl_mat_end_itt)
linesep_vec_end_itt <- rep("", nrow(tbl_mat_end_itt))
is_new_var_end_itt <- c(TRUE, !grepl("^ ", var_names_end_itt[-1]) & var_names_end_itt[-1] != "")
is_last_of_group_end_itt <- c(is_new_var_end_itt[-1], TRUE)
linesep_vec_end_itt[is_last_of_group_end_itt] <- "\\addlinespace"

kbl_out_end_itt <- kable(tbl_mat_end_itt, format = "latex", booktabs = TRUE, escape = TRUE, linesep = linesep_vec_end_itt,
                       caption = "Covariate balance among individuals event-free and uncensored through week 52 (ITT)") |>
  kable_styling(latex_options = c("HOLD_position", "scale_down"))

kbl_str_end_itt <- as.character(kbl_out_end_itt)
kbl_str_end_itt <- gsub("<0.001", "$<$0.001", kbl_str_end_itt, fixed = TRUE)

cat(kbl_str_end_itt)


# pp


tbl_end_pp <- CreateTableOne(
  vars        = vars,
  factorVars  = factor_vars,
  strata      = "treatment_text",
  data        = data_endfu_pp,
  addOverall  = FALSE
)

tbl_mat_end_pp <- print(tbl_end_pp, smd = TRUE, showAllLevels = TRUE, printToggle = FALSE, quote = FALSE, noSpaces = TRUE)
tbl_mat_end_pp <- tbl_mat_end_pp[, colnames(tbl_mat_end_pp) != "test"]

rownames(tbl_mat_end_pp) <- gsub("age_years (mean (SD))", "Age (mean (SD))", rownames(tbl_mat_end_pp), fixed = TRUE)
rownames(tbl_mat_end_pp) <- gsub("bmi (mean (SD))", "BMI (mean (SD))", rownames(tbl_mat_end_pp), fixed = TRUE)
rownames(tbl_mat_end_pp) <- gsub("ndi (mean (SD))", "NDI (mean (SD))", rownames(tbl_mat_end_pp), fixed = TRUE)
rownames(tbl_mat_end_pp) <- gsub("tests_count (mean (SD))", "Test count (mean (SD))", rownames(tbl_mat_end_pp), fixed = TRUE)
rownames(tbl_mat_end_pp) <- gsub("last_vax_infect_weeks (mean (SD))", "Weeks last inf/vax (mean (SD))", rownames(tbl_mat_end_pp), fixed = TRUE)
rownames(tbl_mat_end_pp) <- gsub("sex_admin (%)", "Sex (%)", rownames(tbl_mat_end_pp), fixed = TRUE)
rownames(tbl_mat_end_pp) <- gsub("race (%)", "Race (%)", rownames(tbl_mat_end_pp), fixed = TRUE)
rownames(tbl_mat_end_pp) <- gsub("service_region (%)", "Service region (%)", rownames(tbl_mat_end_pp), fixed = TRUE)
rownames(tbl_mat_end_pp) <- gsub("flu_vax (%)", "Prior year flu (%)", rownames(tbl_mat_end_pp), fixed = TRUE)
rownames(tbl_mat_end_pp) <- gsub("prior_inf (%)", "Prior COVID-19 infection (%)", rownames(tbl_mat_end_pp), fixed = TRUE)
rownames(tbl_mat_end_pp) <- gsub("charlson_cat_fac (%)", "Charlson (%)", rownames(tbl_mat_end_pp), fixed = TRUE)

var_names_end_pp <- rownames(tbl_mat_end_pp)
linesep_vec_end_pp <- rep("", nrow(tbl_mat_end_pp))
is_new_var_end_pp <- c(TRUE, !grepl("^ ", var_names_end_pp[-1]) & var_names_end_pp[-1] != "")
is_last_of_group_end_pp <- c(is_new_var_end_pp[-1], TRUE)
linesep_vec_end_pp[is_last_of_group_end_pp] <- "\\addlinespace"

kbl_out_end_pp <- kable(tbl_mat_end_pp, format = "latex", booktabs = TRUE, escape = TRUE, linesep = linesep_vec_end_pp,
                         caption = "Covariate balance among individuals event-free and uncensored through week 52 (PP)") |>
  kable_styling(latex_options = c("HOLD_position", "scale_down"))

kbl_str_end_pp <- as.character(kbl_out_end_pp)
kbl_str_end_pp <- gsub("<0.001", "$<$0.001", kbl_str_end_pp, fixed = TRUE)

cat(kbl_str_end_pp)


# # Covariate balance by match week -----------------------------------------
# 
# 
# for (wk in weeks) {
#   cat("\n========================================\n")
#   cat(sprintf(" Enrollment Week %d\n", wk))
#   cat("========================================\n")
# 
#   df_wk <- data_Y3 |> filter(index_time == wk)
# 
#   tbl <- CreateTableOne(
#     vars        = vars,
#     factorVars  = factor_vars,
#     strata      = "treatment_text",
#     data        = df_wk,
#     addOverall  = FALSE
#   )
# 
#   tbl_mat <- print(tbl, smd = TRUE, showAllLevels = TRUE, printToggle = FALSE, quote = FALSE, noSpaces = TRUE)
#   tbl_mat <- tbl_mat[, colnames(tbl_mat) != "test"]
# 
#   rownames(tbl_mat) <- gsub("age_years (mean (SD))", "Age (mean (SD))", rownames(tbl_mat), fixed = TRUE)
#   rownames(tbl_mat) <- gsub("bmi (mean (SD))", "BMI (mean (SD))", rownames(tbl_mat), fixed = TRUE)
#   rownames(tbl_mat) <- gsub("ndi (mean (SD))", "NDI (mean (SD))", rownames(tbl_mat), fixed = TRUE)
#   rownames(tbl_mat) <- gsub("tests_count (mean (SD))", "Test count (mean (SD))", rownames(tbl_mat), fixed = TRUE)
#   rownames(tbl_mat) <- gsub("last_vax_infect_weeks (mean (SD))", "Weeks last inf/vax (mean (SD))", rownames(tbl_mat), fixed = TRUE)
#   rownames(tbl_mat) <- gsub("sex_admin (%)", "Sex (%)", rownames(tbl_mat), fixed = TRUE)
#   rownames(tbl_mat) <- gsub("race (%)", "Race (%)", rownames(tbl_mat), fixed = TRUE)
#   rownames(tbl_mat) <- gsub("service_region (%)", "Service region (%)", rownames(tbl_mat), fixed = TRUE)
#   rownames(tbl_mat) <- gsub("flu_vax (%)", "Prior year flu (%)", rownames(tbl_mat), fixed = TRUE)
#   rownames(tbl_mat) <- gsub("prior_inf (%)", "Prior COVID-19 infection (%)", rownames(tbl_mat), fixed = TRUE)
#   rownames(tbl_mat) <- gsub("charlson_cat_fac (%)", "Charlson (%)", rownames(tbl_mat), fixed = TRUE)
# 
#   var_names <- rownames(tbl_mat)
#   linesep_vec <- rep("", nrow(tbl_mat))
# 
#   is_new_var <- c(TRUE, !grepl("^ ", var_names[-1]) & var_names[-1] != "")
#   is_last_of_group <- c(is_new_var[-1], TRUE)
#   linesep_vec[is_last_of_group] <- "\\addlinespace"
# 
#   kbl_out <- kable(tbl_mat, format = "latex", booktabs = TRUE, escape = TRUE, linesep = linesep_vec,
#                    caption = sprintf("Covariate balance, enrollment week %d", wk)) |>
#     kable_styling(latex_options = c("hold_position", "scale_down"))
# 
#   kbl_str <- as.character(kbl_out)
#   kbl_str <- gsub("<0.001", "$<$0.001", kbl_str, fixed = TRUE)
# 
#   cat(kbl_str)
#   cat("\n\n")
# }



