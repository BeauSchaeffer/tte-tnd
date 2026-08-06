##----- Kaiser Causal TTE-TND
##----- Results Table
##----- Intention to Treat
##----- Last updated 2026-08-06


# Packages ----------------------------------------------------------------


library(tidyverse)
library(kableExtra)


# Data --------------------------------------------------------------------


res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results_itt.4/"

# STD Cox
std.itt.cox.pointest <- readRDS(paste0(res_path,"std.itt.cox.pointest.rds"))

# STD Pooled
std.itt.risks.ci <- readRDS(paste0(res_path, "std.itt.risks.ci.rds"))

# TND
tnd.pointest <- readRDS(paste0(res_path,"tnd.itt.pointest.rds"))

# EQC Cox
eqc.itt.HRs.ci <- readRDS(paste0(res_path, "eqc.itt.HRs.ci.rds"))

# EQC Pooled
eqc.itt.risks.ci <- readRDS(paste0(res_path, "eqc.itt.risks.ci.rds"))

# PCI Cox
pci.itt.HRs.ci <- readRDS(paste0(res_path, "pci.itt.HRs.ci.rds"))

# PCI Pooled
pci.itt.risks.ci <- readRDS(paste0(res_path, "pci.itt.risks.ci.rds"))


# Cox estimates table -----------------------------------------------------

label_exposure <- function(term) {
  case_when(
    term == "treatment" ~ "Treatment",
    term == "flu_vax"   ~ "Flu vax",
    TRUE                ~ term
  )
}

std_tbl <- std.itt.cox.pointest |> 
  filter(term %in% c("treatment", "flu_vax")) |> 
  transmute(
    approach  = "STD",
    term,
    exposure  = label_exposure(term),
    effect    = estimate,
    ci_low    = conf.low,
    ci_high   = conf.high
  ) |> 
  select(-term)

tnd_tbl <- tnd.pointest |> 
  filter(term %in% c("treatment", "flu_vax")) |> 
  transmute(
    approach = "TND",
    term,
    exposure = label_exposure(term),
    effect   = OR,
    ci_low   = LCL,
    ci_high  = UCL
  ) |> 
  select(-term)


eqc_tbl <- eqc.itt.HRs.ci %>%
  pivot_longer(everything(), names_to = "name", values_to = "value") %>%
  mutate(
    exposure = case_when(
      grepl("^treatHR",  name) ~ "Treatment",
      grepl("^fluvaxHR", name) ~ "Flu vax",
      TRUE ~ NA_character_
    ),
    metric = case_when(
      grepl("_lo$", name) ~ "ci_low",
      grepl("_hi$", name) ~ "ci_high",
      name %in% c("treatHR", "fluvaxHR") ~ "effect",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(exposure), !is.na(metric)) |> 
  select(exposure, metric, value) |> 
  pivot_wider(names_from = metric, values_from = value) |> 
  mutate(approach = "EQC") |> 
  select(approach, exposure, effect, ci_low, ci_high)

pci_tbl <- pci.itt.HRs.ci |> 
  transmute(
    approach = "PCI",
    exposure = "Treatment",
    effect   = treatHR,
    ci_low   = treatHR_lo,
    ci_high  = treatHR_hi
  )

draft_table <- bind_rows(std_tbl, tnd_tbl, eqc_tbl, pci_tbl) |> 
  mutate(ve = (1-effect)*100,
         ve_low = (1-ci_high)*100,
         ve_high = (1-ci_low)*100) |> 
  mutate(
    effect.ci = sprintf("%.2f (%.2f, %.2f)", effect, ci_low, ci_high),
    ve.ci = sprintf("%.2f (%.2f, %.2f)", ve, ve_low, ve_high)
  ) |> 
  select(approach, exposure, effect.ci, ve.ci)

draft_table

wide_table <- draft_table |>
  
  pivot_longer(
    cols = c(effect.ci, ve.ci),
    names_to = "metric",
    values_to = "value"
  ) |>

  mutate(
    metric = recode(metric,
                    "effect.ci" = "HR",
                    "ve.ci"     = "VE (%)"
    ),
    metric = factor(metric, levels = c("HR", "VE (%)"))
  ) |>
  
  pivot_wider(
    names_from  = c(exposure, metric),
    values_from = value,
    names_sep   = "_"
  ) |>
  
  select(
    approach,
    Treatment_HR, `Treatment_VE (%)`,
    `Flu vax_HR`, `Flu vax_VE (%)`
  )

wide_table


wide_table |>
  kbl(
    format = "latex",
    caption   = "Cox results",
    label     = "cox_res_itt",
    booktabs = TRUE,
    align = "lcccc",
    col.names = c("Approach", "HR", "VE (%)", "HR", "VE (%)")
  ) |>
  add_header_above(c(" " = 1, "Bivalent Booster" = 2, "Prior year flu vaccine" = 2)) |>
  kable_styling(full_width = FALSE)


# Pooled estimates table --------------------------------------------------

# horizons <- c(1, 8, 24, 40, 52)
horizons <- seq(1,52,1)

make_pooled_risk_table <- function(df, approach, horizons) {
  df |> 
    filter(time_end %in% horizons) %>%
    mutate(
      `Risk Ratio (95\\% CI)` = sprintf("%.2f (%.2f, %.2f)", rr, rr_lo, rr_hi),
      `Risk Difference (95\\% CI)` = sprintf("%.4f (%.4f, %.4f)", rd, rd_lo, rd_hi),
    ) |> 
    transmute(
      `Time (weeks)` = time_end,
      `Risk Ratio (95\\% CI)`,
      `Risk Difference (95\\% CI)`
      )
}

std_risk_tbl <- make_pooled_risk_table(std.itt.risks.ci, "STD", horizons)
eqc_risk_tbl <- make_pooled_risk_table(eqc.itt.risks.ci, "EQC", horizons)
pci_risk_tbl <- make_pooled_risk_table(pci.itt.risks.ci, "PCI", horizons)

knitr::kable(
  std_risk_tbl,
  format = "latex",
  booktabs = TRUE,
  caption = "Pooled logistic results by week (STD)",
  label = "std_rr_rd_itt",
  escape = F
)

knitr::kable(
  eqc_risk_tbl,
  format = "latex",
  booktabs = TRUE,
  caption = "Pooled logistic results by week (EQC)",
  label = "eqc_rr_rd_itt",
  escape = F
)

knitr::kable(
  pci_risk_tbl,
  format = "latex",
  booktabs = TRUE,
  caption = "Pooled logistic results by week (PCI)",
  label = "pci_rr_rd_itt",
  escape = F
)










