#!/bin/bash
#SBATCH -J tnd_itt
#SBATCH -o logs/tnd_itt_%j.out
#SBATCH -e logs/tnd_itt_%j.err
#SBATCH -p shared
#SBATCH -c 1
#SBATCH -t 0-00:10
#SBATCH --mem=10000 # in Mb

# RStudio Server container + package library
R_LIB="/n/home07/bschaeffer/R/ifxrstudio/RELEASE_3_20"
RSTUDIO_IMG="/n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_20.sif"

singularity exec \
	--cleanenv \
	--env R_LIBS_USER=${R_LIB} \
	${RSTUDIO_IMG} \
	Rscript scripts/itt_scripts/tnd_analysis_itt.R

# run with
# sbatch slurm/tnd_itt.sbatch
