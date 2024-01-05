#!/bin/bash

#SBATCH --time=12:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=80   # 50 processor core(s) per node X 2 threads per core
#SBATCH --mem=1200G   # maximum memory per node
#SBATCH --partition=priority-mem    # large-memory node(s)
#SBATCH --qos msn
#SBATCH --job-name="polyRAD"
#SBATCH --mail-user=william.rolling@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output="polyRADout" # job standard output file (%j replaced by job id)
#SBATCH --error="polyRAD.err" # job standard error file (%j replaced by job id)

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load r/4.2.3
R CMD BATCH --vanilla polyRAD.R
