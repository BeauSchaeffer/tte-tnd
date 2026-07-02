#!/bin/bash
#SBATCH -J std_pp_boot
#SBATCH -o logs/std_pp_boot_%j.out
#SBATCH -e logs/std_pp_boot_%j.err
#SBATCH -p shared
#SBATCH -c 1
#SBATCH -t 0-12:00
#SBATCH --mem=64000 # in Mb

# RStudio Server container + package library
R_LIB="/n/home07/bschaeffer/R/ifxrstudio/RELEASE_3_20"
RSTUDIO_IMG="/n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_20.sif"

singularity exec \
	--cleanenv \
	--env R_LIBS_USER=${R_LIB} \
	${RSTUDIO_IMG} \
	Rscript scripts/per_prot_nocw/standard_analysis_pooled_pp_boot.R

# run with
# sbatch slurm/std_pp_boot.sbatch
