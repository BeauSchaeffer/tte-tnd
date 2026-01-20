##----- Kaiser Causal TTE-TND
##----- Results Table


# Packages ----------------------------------------------------------------


library(tidyverse)


# Data --------------------------------------------------------------------


res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results/"

# STD Cox
std.itt.cox.pointest <- readRDS(paste0(res_path,"std.itt.cox.pointest.rds"))

# STD Pooled
std.itt.risks.ci <- readRDS(paste0(res_path, "std.itt.risks.ci.rds"))

# TND
tnd.pointest <- readRDS(paste0(res_path,"tnd.pointest.rds"))

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

draft_table <- bind_rows(std_tbl, tnd_tbl, eqc_tbl, pci_tbl) %>%
  mutate(
    effect.ci = sprintf("%.2f (%.2f–%.2f)", effect, ci_low, ci_high)
  ) %>%
  select(approach, exposure, effect.ci)

draft_table

wide_table <- draft_table |> 
  mutate(
    exposure = case_when(
      exposure %in% c("Flu vax", "Flu vaccine") ~ "Flu vaccine",
      TRUE ~ exposure
    )
  ) |> 
  pivot_wider(
    names_from  = exposure,
    values_from = effect.ci
  )

wide_table

knitr::kable(
  wide_table,
  col.names = c(
    "Approach",
    "Treatment (95% CI)",
    "Flu vaccine (95% CI)"
  ),
  align = c("l", "c", "c"))

knitr::kable(
  wide_table,
  format = "latex",
  booktabs = TRUE,
  col.names = c(
    "Approach",
    "Treatment (95\\% CI)",
    "Flu vaccine (95\\% CI)"
  ),
  align = c("l", "c", "c"))


# Pooled estimates table --------------------------------------------------

horizons <- c(4, 24, 52)

make_risk_summary <- function(df, approach, horizons) {
  df %>%
    filter(time_end %in% horizons) %>%
    transmute(
      approach = approach,
      time_end,
      `Risk (no booster)` = sprintf("%.4f (%.4f–%.4f)", risk0, risk0_lo, risk0_hi),
      `Risk (booster)` = sprintf("%.4f (%.4f–%.4f)", risk1, risk1_lo, risk1_hi),
      `Risk Ratio (booster / no booster)` = sprintf(
        "%.4f", risk1/risk0
      )
    )
}

risk_long <- bind_rows(
  make_risk_summary(std.itt.risks.ci, "STD", horizons),
  make_risk_summary(eqc.itt.risks.ci, "EQC", horizons),
  make_risk_summary(pci.itt.risks.ci, "PCI", horizons)
)

risk_wide <- risk_long %>%
  pivot_longer(
    cols = c(`Risk (no booster)`, `Risk (booster)`, `Risk Ratio (booster / no booster)`),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(col = paste0(metric, " @ t=", time_end)) %>%
  select(approach, col, value) %>%
  pivot_wider(names_from = col, values_from = value) %>%
  arrange(approach)

knitr::kable(
  risk_wide,
  format = "latex",
  booktabs = TRUE,
  caption = "Pooled logistic estimated cumulative risks by approach and time horizon",
  align = c("l", rep("c", ncol(risk_wide) - 1))
)

risk_wide2 <- risk_long %>%
  rename(
    `No Booster` = `Risk (no booster)`,
    `Booster` = `Risk (booster)`,
    RR = `Risk Ratio (booster / no booster)`
  ) %>%
  pivot_longer(cols = c(`No Booster`, Booster, RR),
               names_to = "metric", values_to = "value") %>%
  mutate(col = paste0(metric, " t=", time_end)) %>%
  select(approach, col, value) %>%
  pivot_wider(names_from = col, values_from = value) %>%
  arrange(approach)

knitr::kable(
  risk_wide2,
  format = "latex",
  booktabs = TRUE,
  caption = "Pooled logistic cumulative risk (95\\% CI) and risk difference by approach and time horizon",
  align = c("l", rep("c", ncol(risk_wide2) - 1))
)

knitr::kable(
  risk_wide2,
  align = c("l", rep("c", ncol(risk_wide2) - 1))
)


