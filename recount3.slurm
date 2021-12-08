#!/bin/bash

#sbatch --job-name=recount3
#sbatch --mem=256G --cpus-per-task=32
#sbatch --time=5:00:00:00
#sbatch --qos=1week
#SBATCH --output=out/recount3.out
#SBATCH --error=out/recount3.err

module load R/4.1.1
module load gsl/2.5
module load gcc/7.2.0
cd /Genomics/ayroleslab2/diogro/projects/exp_var
Rscript runRecount.R
