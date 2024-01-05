#!/bin/bash


#SBATCH --time=04:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=96   # 48 processor core(s) per node X 2 threads per core
#SBATCH --mem=6400M   # maximum memory per node
#SBATCH --qos msn
#SBATCH --partition=priority    # standard node(s)
#SBATCH --job-name="Admixture"
#SBATCH --output="Admix.out" # job standard output file (%j replaced by job id)
#SBATCH --error="Admix.error" # job standard error file (%j replaced by job id)

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load admixture/1.3.0
admixture ../plinkGenos/GenotypeFile.bed 1
admixture ../plinkGenos/GenotypeFile.bed 2
admixture ../plinkGenos/GenotypeFile.bed 3
admixture ../plinkGenos/GenotypeFile.bed 4
admixture ../plinkGenos/GenotypeFile.bed 5
admixture ../plinkGenos/GenotypeFile.bed 6
admixture ../plinkGenos/GenotypeFile.bed 7
admixture ../plinkGenos/GenotypeFile.bed 8
admixture ../plinkGenos/GenotypeFile.bed 9
admixture ../plinkGenos/GenotypeFile.bed 10
