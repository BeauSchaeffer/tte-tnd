#!/bin/bash
#SBATCH -J pci_pp_cox_with_boot
#SBATCH -o logs/pci_pp_cox_with_boot_%j.out
#SBATCH -e logs/pci_pp_cox_with_boot_%j.err
#SBATCH -p shared
#SBATCH -c 1
#SBATCH -t 0-1:00
#SBATCH --mem=64000 # in Mb

# RStudio Server container + package library
R_LIB="/n/home07/bschaeffer/R/ifxrstudio/RELEASE_3_20"
RSTUDIO_IMG="/n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_20.sif"

singularity exec \
	--cleanenv \
	--env R_LIBS_USER=${R_LIB} \
	${RSTUDIO_IMG} \
	Rscript scripts/per_prot_nocw/proximal_inference_analysis_pp.R

# run with
# sbatch slurm/pci_pp_cox_with_boot.sbatch 
