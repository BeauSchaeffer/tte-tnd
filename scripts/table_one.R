##----- Beau Schaeffer
##----- Kaiser TTE-TND
##----- Table 1(s)

library(tidyverse)
library(tableone)
library(kableExtra)

# Load data ---------------------------------------------------------------

data_Y3 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch/data_Y3_weekmatch.rds")

# Table 1 by Week ---------------------------------------------------------



vars <- c("age_years", "bmi", "ndi", "tests_count", "last_vax_infect_weeks",
          "sex_admin", "race", "service_region", "flu_vax",
          "prior_inf", "charlson_cat_fac")

factor_vars <- c("sex_admin", "race", "service_region", "flu_vax",
                 "prior_inf", "charlson_cat_fac")

weeks <- sort(unique(data_Y3$index_time))

data_Y3 <- data_Y3 |>
  mutate(treatment_text = factor(treatment, levels = c(1, 0),
                                 labels = c("Booster", "No Booster")))

for (wk in weeks) {
  cat("\n========================================\n")
  cat(sprintf(" Enrollment Week %d\n", wk))
  cat("========================================\n")
  
  df_wk <- data_Y3 |> filter(index_time == wk)
  
  tbl <- CreateTableOne(
    vars        = vars,
    factorVars  = factor_vars,
    strata      = "treatment_text",
    data        = df_wk,
    addOverall  = FALSE
  )
  
  tbl_mat <- print(tbl, smd = TRUE, showAllLevels = TRUE, printToggle = FALSE, quote = FALSE, noSpaces = TRUE)
  tbl_mat <- tbl_mat[, colnames(tbl_mat) != "test"]
  
  rownames(tbl_mat) <- gsub("age_years (mean (SD))", "Age (mean (SD))", rownames(tbl_mat), fixed = TRUE)
  rownames(tbl_mat) <- gsub("bmi (mean (SD))", "BMI (mean (SD))", rownames(tbl_mat), fixed = TRUE)
  rownames(tbl_mat) <- gsub("ndi (mean (SD))", "NDI (mean (SD))", rownames(tbl_mat), fixed = TRUE)
  rownames(tbl_mat) <- gsub("tests_count (mean (SD))", "Test count (mean (SD))", rownames(tbl_mat), fixed = TRUE)
  rownames(tbl_mat) <- gsub("last_vax_infect_weeks (mean (SD))", "Weeks last inf/vax (mean (SD))", rownames(tbl_mat), fixed = TRUE)
  rownames(tbl_mat) <- gsub("sex_admin (%)", "Sex (%)", rownames(tbl_mat), fixed = TRUE)
  rownames(tbl_mat) <- gsub("race (%)", "Race (%)", rownames(tbl_mat), fixed = TRUE)
  rownames(tbl_mat) <- gsub("service_region (%)", "Service region (%)", rownames(tbl_mat), fixed = TRUE)
  rownames(tbl_mat) <- gsub("flu_vax (%)", "Prior year flu (%)", rownames(tbl_mat), fixed = TRUE)
  rownames(tbl_mat) <- gsub("prior_inf (%)", "Prior COVID-19 infection (%)", rownames(tbl_mat), fixed = TRUE)
  rownames(tbl_mat) <- gsub("charlson_cat_fac (%)", "Charlson (%)", rownames(tbl_mat), fixed = TRUE)
  
  var_names <- rownames(tbl_mat)
  linesep_vec <- rep("", nrow(tbl_mat))
  
  is_new_var <- c(TRUE, !grepl("^ ", var_names[-1]) & var_names[-1] != "")
  is_last_of_group <- c(is_new_var[-1], TRUE)
  linesep_vec[is_last_of_group] <- "\\addlinespace"
  
  kbl_out <- kable(tbl_mat, format = "latex", booktabs = TRUE, escape = TRUE, linesep = linesep_vec,
                   caption = sprintf("Covariate balance, enrollment week %d", wk)) |>
    kable_styling(latex_options = c("hold_position", "scale_down"))
  
  kbl_str <- as.character(kbl_out)
  kbl_str <- gsub("<0.001", "$<$0.001", kbl_str, fixed = TRUE)
  
  cat(kbl_str)
  cat("\n\n")
}


# Overall Table 1 ---------------------------------------------------------

# Overall Table 1 ---------------------------------------------------------

tbl_overall <- CreateTableOne(
  vars        = vars,
  factorVars  = factor_vars,
  strata      = "treatment_text",
  data        = data_Y3,
  addOverall  = FALSE
)

tbl_mat_overall <- print(tbl_overall, smd = TRUE, showAllLevels = TRUE, printToggle = FALSE, quote = FALSE, noSpaces = TRUE)
tbl_mat_overall <- tbl_mat_overall[, colnames(tbl_mat_overall) != "test"]

rownames(tbl_mat_overall) <- gsub("age_years (mean (SD))", "Age (mean (SD))", rownames(tbl_mat_overall), fixed = TRUE)
rownames(tbl_mat_overall) <- gsub("bmi (mean (SD))", "BMI (mean (SD))", rownames(tbl_mat_overall), fixed = TRUE)
rownames(tbl_mat_overall) <- gsub("ndi (mean (SD))", "NDI (mean (SD))", rownames(tbl_mat_overall), fixed = TRUE)
rownames(tbl_mat_overall) <- gsub("tests_count (mean (SD))", "Test count (mean (SD))", rownames(tbl_mat_overall), fixed = TRUE)
rownames(tbl_mat_overall) <- gsub("last_vax_infect_weeks (mean (SD))", "Weeks last inf/vax (mean (SD))", rownames(tbl_mat_overall), fixed = TRUE)
rownames(tbl_mat_overall) <- gsub("sex_admin (%)", "Sex (%)", rownames(tbl_mat_overall), fixed = TRUE)
rownames(tbl_mat_overall) <- gsub("race (%)", "Race (%)", rownames(tbl_mat_overall), fixed = TRUE)
rownames(tbl_mat_overall) <- gsub("service_region (%)", "Service region (%)", rownames(tbl_mat_overall), fixed = TRUE)
rownames(tbl_mat_overall) <- gsub("flu_vax (%)", "Prior year flu (%)", rownames(tbl_mat_overall), fixed = TRUE)
rownames(tbl_mat_overall) <- gsub("prior_inf (%)", "Prior COVID-19 infection (%)", rownames(tbl_mat_overall), fixed = TRUE)
rownames(tbl_mat_overall) <- gsub("charlson_cat_fac (%)", "Charlson (%)", rownames(tbl_mat_overall), fixed = TRUE)

var_names_overall <- rownames(tbl_mat_overall)
linesep_vec_overall <- rep("", nrow(tbl_mat_overall))
is_new_var_overall <- c(TRUE, !grepl("^ ", var_names_overall[-1]) & var_names_overall[-1] != "")
is_last_of_group_overall <- c(is_new_var_overall[-1], TRUE)
linesep_vec_overall[is_last_of_group_overall] <- "\\addlinespace"

kbl_out_overall <- kable(tbl_mat_overall, format = "latex", booktabs = TRUE, escape = TRUE, linesep = linesep_vec_overall,
                         caption = "Covariate balance, all weeks") |>
  kable_styling(latex_options = c("hold_position", "scale_down"))

kbl_str_overall <- as.character(kbl_out_overall)
kbl_str_overall <- gsub("<0.001", "$<$0.001", kbl_str_overall, fixed = TRUE)

cat(kbl_str_overall)
