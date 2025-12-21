#!/bin/bash
#SBATCH --job-name=std_itt_boot
#SBATCH --partition=shared
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --output=/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/logs/std_itt_boot_%j.out
#SBATCH --error=/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/logs/std_itt_boot_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=YOUR_EMAIL@hsph.harvard.edu

set -euo pipefail

# ---- (optional) make sure a logs dir exists ----
mkdir -p /n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/logs

# ---- load R (edit to match your cluster) ----
module purge
module load R/4.3.2

# ---- run ----
Rscript --vanilla /n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/scripts/std_itt_bootstrap.R

