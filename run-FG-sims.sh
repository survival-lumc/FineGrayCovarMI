#!/bin/bash
#SBATCH -J MI-FG
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --partition=all
#SBATCH --time=0-00:59:00
#SBATCH --mem-per-cpu=5G


# Purge
module purge

# Load module
module load statistical/R/4.3.1/gcc.8.5.0

# Run imps
Rscript -e 'targets::tar_make_future(workers = parallelly::availableCores(logical = FALSE))'
