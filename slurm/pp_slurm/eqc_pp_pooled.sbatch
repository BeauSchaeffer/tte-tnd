#!/bin/bash
#SBATCH -J eqc_pp_pooled
#SBATCH -o logs/eqc_pp_pooled_%j.out
#SBATCH -e logs/eqc_pp_pooled_%j.err
#SBATCH -p shared
#SBATCH -c 1
#SBATCH -t 0-0:30
#SBATCH --mem=64000 # in Mb

# RStudio Server container + package library
R_LIB="/n/home07/bschaeffer/R/ifxrstudio/RELEASE_3_20"
RSTUDIO_IMG="/n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_20.sif"

singularity exec \
	--cleanenv \
	--env R_LIBS_USER=${R_LIB} \
	${RSTUDIO_IMG} \
	Rscript scripts/pp_scripts/equi_confounding_pooled_pp.R

# run with
# sbatch slurm/pp_slurm/eqc_pp_pooled.sbatch
