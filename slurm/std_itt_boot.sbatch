#!/bin/bash
#SBATCH -J std_itt_boot
#SBATCH -o logs/std_itt_boot_%j.out
#SBATCH -e logs/std_itt_boot_%j.err
#SBATCH -p shared
#SBATCH -c 1
#SBATCH -t 0-01:00
#SBATCH --mem=64000 # in Mb

# Load required software modules
module load R

# Run program
R CMD BATCH --no-save --no-restore \
  scripts/standard_analysis_pooled_itt_boot.R \
  logs/standard_analysis_pooled_itt_boot.Rout
