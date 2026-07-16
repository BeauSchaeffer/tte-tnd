#!/bin/bash
#SBATCH -J pci_itt_boot
#SBATCH -o logs/pci_itt_boot_%A_%a.out
#SBATCH -e logs/pci_itt_boot_%A_%a.err
#SBATCH -p shared
#SBATCH -c 1
#SBATCH -t 0-45:00
#SBATCH --mem=120000
#SBATCH --array=1-50

# RStudio Server container + package library
R_LIB="/n/home07/bschaeffer/R/ifxrstudio/RELEASE_3_20"
RSTUDIO_IMG="/n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_20.sif"

# SLURM_ARRAY_TASK_ID is 1..50

singularity exec \
  --cleanenv \
  --env R_LIBS_USER=${R_LIB} \
  ${RSTUDIO_IMG} \
  Rscript scripts/itt_scripts/proximal_inference_pooled_itt_boot_onerep.R ${SLURM_ARRAY_TASK_ID}

# run with
# sbatch slurm/itt_slurm/pci_itt_pooled_boot.sbatch