#!/bin/bash
#SBATCH -J netB_pci_pp_pooled
#SBATCH -o logs/netB_pci_pp_pooled_%j.out
#SBATCH -e logs/netB_pci_pp_pooled_%j.err
#SBATCH -p shared
#SBATCH -c 1
#SBATCH -t 0-01:30
#SBATCH --mem=160000 # in Mb

# RStudio Server container + package library
R_LIB="/n/home07/bschaeffer/R/ifxrstudio/RELEASE_3_20"
RSTUDIO_IMG="/n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_20.sif"

singularity exec \
	--cleanenv \
	--env R_LIBS_USER=${R_LIB} \
	${RSTUDIO_IMG} \
	Rscript scripts/net_risk/proximal_inference_pooled_netB_pp.R

# run with
# sbatch slurm/net_risk_slurm/netB_pci_pp_pooled.sbatch
