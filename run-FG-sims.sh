#!/bin/bash
#SBATCH -J MI-FG
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=short
#SBATCH --time=0-00:30:00 
#SBATCH --mem=5G


# Purge
module purge

# Load module
module load statistical/R/4.2.1

# Run imps
Rscript -e 'targets::tar_make()' #./one-replication_cens.R
