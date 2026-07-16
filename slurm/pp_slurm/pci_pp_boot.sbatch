#!/bin/bash
#SBATCH -J pci_pp_boot
#SBATCH -o logs/pci_pp_boot_%A_%a.out
#SBATCH -e logs/pci_pp_boot_%A_%a.err
#SBATCH -p shared
#SBATCH -c 1
#SBATCH -t 0-30:00
#SBATCH --mem=100000
#SBATCH --array=1-100

# RStudio Server container + package library
R_LIB="/n/home07/bschaeffer/R/ifxrstudio/RELEASE_3_20"
RSTUDIO_IMG="/n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_20.sif"

# SLURM_ARRAY_TASK_ID is 1..100

singularity exec \
  --cleanenv \
  --env R_LIBS_USER=${R_LIB} \
  ${RSTUDIO_IMG} \
  Rscript scripts/per_prot_nocw/proximal_inference_pooled_pp_boot_onerep.R ${SLURM_ARRAY_TASK_ID}

# run with
# sbatch slurm/pci_pp_boot.sbatch