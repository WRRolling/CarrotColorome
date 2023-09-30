#!/bin/bash

#SBATCH --time=04:0:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=80   # 36 processor core(s) per node X 2 threads per core
#SBATCH --mem=1494G   # maximum memory per node
#SBATCH --partition=scavenger    # standard node(s)
#SBATCH --job-name="Filter.2"
#SBATCH --mail-user=william.rolling@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output="Filter.2" # job standard output file (%j replaced by job id)
#SBATCH --error="Filter.2" # job standard error file (%j replaced by job id)

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

# load vcftools module
module load vcftools/0.1.16
module load samtools/1.17
# options for how to filter VCF file

perl ../../bb.vcf --infile="FiltWI18.vcf.gz" --outfile="FiltWI18.2.vcf.gz" --task=het --minratio=0.3 --maxratio=0.7
