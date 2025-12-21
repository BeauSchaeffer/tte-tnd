##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Standard Analysis Pooled


# Packages ----------------------------------------------------------------


library(tidyverse)
library(tableone)
library(gridExtra)
library(ragg)


# Data --------------------------------------------------------------------


data_Y2 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_Y2.rds")


# Clean up names ----------------------------------------------------------


data_Y2_clean <- data_Y2 |> 
  rename(Sex=sex_admin,
         Age=age_years,
         NDI=ndi,
         BMI=bmi,
         `Prior year influenza vax`=flu_vax,
         `Prior COVID-19 infection`=prior_inf,
         `Total COVID-19 test count`=tests_count,
         Race=race,
         `Service region`=service_region,
         `Weeks since last vax or infection`=last_vax_infect_weeks,
         `Charlson index`=charlson_cat_fac,
         Treatment=treatment) |> 
  
  mutate(Sex=case_when(Sex=="F" ~ "Female",
                       Sex=="M" ~"Male"),
         
         `Prior year influenza vax`=case_when(
           `Prior year influenza vax`==1~"Yes",
           `Prior year influenza vax`==0~"No"),
         
         `Prior COVID-19 infection`=case_when(
           `Prior COVID-19 infection`==1~"Yes",
           `Prior COVID-19 infection`==0~"No"),
         
         Treatment=case_when(Treatment==1~"Booster",
                             Treatment==0~"No Booster"),
         
         Race=case_when(Race=="Unknown/ot"~"Unknown/other",
                        TRUE~Race)
         )

# Table One ---------------------------------------------------------------

vars <- c(
  "Sex",
  "Age",
  "NDI",
  "BMI",
  "Prior year influenza vax",
  "Prior COVID-19 infection",
  "Total COVID-19 test count",
  "Race",
  "Service region",
  "Weeks since last vax or infection",
  "Charlson index"
)

factor_vars <- c(
  "Sex",
  "Prior year influenza vax",
  "Prior COVID-19 infection",
  "Race",
  "Service region",
  "Charlson index"
)

table1 <- CreateTableOne(
  vars        = vars,
  strata      = "Treatment",
  data        = data_Y2_clean,
  factorVars  = factor_vars,
  includeNA   = TRUE
)

tab_mat <- print(
  table1,
  showAllLevels = TRUE,
  quote = FALSE,
  noSpaces = TRUE,
  printToggle = FALSE
)

agg_png("figures/table1.png", width = 2000, height = 1600, res = 200)
grid.table(tab_mat)
dev.off()

