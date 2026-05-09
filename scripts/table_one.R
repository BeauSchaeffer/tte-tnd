##----- Beau Schaeffer
##----- Kaiser TTE-TND
##----- Table 1(s)

library(tidyverse)
library(tableone)
library(kableExtra)

# Load data ---------------------------------------------------------------

data_Y3 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch/data_Y3_weekmatch.rds")


# Full Table 1 ------------------------------------------------------------


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
  cat(sprintf("  Enrollment Week %d\n", wk))
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
  
  kbl_out <- kable(tbl_mat, format = "latex", booktabs = TRUE,
                   caption = sprintf("Covariate balance, enrollment week %d", wk)) |>
    kable_styling(latex_options = c("hold_position", "scale_down")) |>
    add_footnote("SMD: standardized mean difference", notation = "none")
  
  cat(as.character(kbl_out))
  cat("\n\n")
}


# Full table 1 ------------------------------------------------------------

tbl <- CreateTableOne(
  vars        = vars,
  factorVars  = factor_vars,
  strata      = "treatment_text",
  data        = data_Y3,
  addOverall  = FALSE
)

tbl_mat <- print(tbl, smd = TRUE, showAllLevels = TRUE, printToggle = FALSE, quote = FALSE, noSpaces = TRUE)
tbl_mat <- apply(tbl_mat, 2, function(x) gsub("<", "$<$", x, fixed = TRUE))
tbl_mat <- apply(tbl_mat, 2, function(x) gsub(">", "$>$", x, fixed = TRUE))
tbl_mat <- apply(tbl_mat, 2, function(x) gsub("%", "\\%", x, fixed = TRUE))

kbl_out <- kable(tbl_mat, format = "latex", booktabs = TRUE, escape=FALSE,
                 caption = sprintf("Covariate balance, all enrollment weeks", wk)) |>
  kable_styling(latex_options = c("hold_position", "scale_down")) |>
  add_footnote("SMD: standardized mean difference", notation = "none")

cat(as.character(kbl_out))
