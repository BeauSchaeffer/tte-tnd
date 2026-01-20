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

head(std.itt.cox.pointest)

head(tnd.pointest)

head(eqc.itt.HRs.ci)

head(pci.itt.HRs.ci)





