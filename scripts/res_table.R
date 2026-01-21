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

horizons <- c(2, 8, 16, 24, 52)

make_pooled_risk_table <- function(df, approach, horizons) {
  df %>%
    filter(time_end %in% horizons) %>%
    mutate(
      `Risk A=0 (95\\% CI)` = sprintf("%.4f (%.4f--%.4f)", risk0, risk0_lo, risk0_hi),
      `Risk A=1 (95\\% CI)` = sprintf("%.4f (%.4f--%.4f)", risk1, risk1_lo, risk1_hi),
      RR = risk1 / risk0,
      RD = risk1 - risk0
    ) %>%
    transmute(
      `t (week)` = time_end,
      `Risk A=0 (95\\% CI)`,
      `Risk A=1 (95\\% CI)`,
      RR = sprintf("%.3f", RR),
      RD = sprintf("%.4f", RD)
    )
}

std_risk_tbl <- make_pooled_risk_table(std.itt.risks.ci, "STD", horizons)
eqc_risk_tbl <- make_pooled_risk_table(eqc.itt.risks.ci, "EQC", horizons)
pci_risk_tbl <- make_pooled_risk_table(pci.itt.risks.ci, "PCI", horizons)

knitr::kable(
  std_risk_tbl,
  format = "latex",
  booktabs = TRUE,
  caption = "Pooled logistic cumulative risks by time horizon (STD)"
)

knitr::kable(
  eqc_risk_tbl,
  format = "latex",
  booktabs = TRUE,
  caption = "Pooled logistic cumulative risks by time horizon (EQC)"
)

knitr::kable(
  pci_risk_tbl,
  format = "latex",
  booktabs = TRUE,
  caption = "Pooled logistic cumulative risks by time horizon (PCI)"
)










