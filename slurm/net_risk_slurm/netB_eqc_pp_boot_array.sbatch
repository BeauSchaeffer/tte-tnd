#!/bin/bash
#SBATCH -J netB_eqc_pp_boot
#SBATCH -o logs/netB_eqc_pp_boot_logs/netB_eqc_pp_boot_%A_%a.out
#SBATCH -e logs/netB_eqc_pp_boot_logs/netB_eqc_pp_boot_%A_%a.err
#SBATCH -p shared
#SBATCH -c 1
#SBATCH -t 0-00:30
#SBATCH --mem=64000
#SBATCH --array=1-100

# RStudio Server container + package library
R_LIB="/n/home07/bschaeffer/R/ifxrstudio/RELEASE_3_20"
RSTUDIO_IMG="/n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_20.sif"

# SLURM_ARRAY_TASK_ID is 1..100

singularity exec \
  --cleanenv \
  --env R_LIBS_USER=${R_LIB} \
  ${RSTUDIO_IMG} \
  Rscript scripts/net_risk/equi_confounding_pooled_netB_pp_boot_onerep.R ${SLURM_ARRAY_TASK_ID}

# run with
# sbatch slurm/net_risk_slurm/netB_eqc_pp_boot_array.sbatch
