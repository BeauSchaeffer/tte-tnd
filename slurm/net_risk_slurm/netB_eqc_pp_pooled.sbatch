#!/bin/bash
#SBATCH -J netB_eqc_pp_pooled
#SBATCH -o logs/netB_eqc_pp_pooled_%j.out
#SBATCH -e logs/netB_eqc_pp_pooled_%j.err
#SBATCH -p shared
#SBATCH -c 1
#SBATCH -t 0-1:30
#SBATCH --mem=64000 # in Mb

# RStudio Server container + package library
R_LIB="/n/home07/bschaeffer/R/ifxrstudio/RELEASE_3_20"
RSTUDIO_IMG="/n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_20.sif"

singularity exec \
	--cleanenv \
	--env R_LIBS_USER=${R_LIB} \
	${RSTUDIO_IMG} \
	Rscript scripts/net_risk/equi_confounding_pooled_netB_pp.R

# run with
# sbatch slurm/net_risk_slurm/netB_eqc_pp_pooled.sbatch
