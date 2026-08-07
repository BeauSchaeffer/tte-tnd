#!/bin/bash
#SBATCH -J pci_pp_pooled_add
#SBATCH -o logs/pci_pp_pooled_add_%j.out
#SBATCH -e logs/pci_pp_pooled_add_%j.err
#SBATCH -p shared
#SBATCH -c 1
#SBATCH -t 0-0:45
#SBATCH --mem=120000 # in Mb

# RStudio Server container + package library
R_LIB="/n/home07/bschaeffer/R/ifxrstudio/RELEASE_3_20"
RSTUDIO_IMG="/n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_20.sif"

singularity exec \
	--cleanenv \
	--env R_LIBS_USER=${R_LIB} \
	${RSTUDIO_IMG} \
	Rscript scripts/pp_scripts/proximal_inference_pooled_pp_additive.R

# run with
# sbatch slurm/pp_slurm/pci_pp_pooled_add.sbatch
