#!/bin/bash

#sbatch --job-name=recount3
#sbatch --mem=450G --cpus-per-task=64
#sbatch --time=5:00:00:00
#sbatch --qos=1week
#SBATCH --output=out/recount3_2.out
#SBATCH --error=out/recount3_2.err

module load R/4.1.1
cd /Genomics/ayroleslab2/diogro/projects/exp_var
Rscript runRecount.R
