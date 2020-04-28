#!/bin/bash
#SBATCH --job-name=drake
#SBATCH --partition=compute
#SBATCH --time=84:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --input=none
#SBATCH --output=drake_%A_%a.out
#SBATCH --error=drake_%A_%a.err

module load R/3.6.1 # Uncomment if R is an environment module.
R --no-save CMD BATCH make.R
