#!/bin/bash
#SBATCH -J MI-FG
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --partition=all
#SBATCH --time=0-12:00:00
#SBATCH --mem-per-cpu=10G


# Purge
module purge

# Load module
module load statistical/R/4.2.1

# Run imps
Rscript -e 'targets::tar_make_future(workers = parallelly::availableCores(logical = FALSE))'
