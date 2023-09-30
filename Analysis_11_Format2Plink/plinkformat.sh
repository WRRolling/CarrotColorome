#!/bin/bash


#SBATCH --time=01:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=40   # 50 processor core(s) per node X 2 threads per core
#SBATCH --qos msn
#SBATCH --partition=priority
#SBATCH --job-name="Create Plink Files"
#SBATCH --mail-user=william.rolling@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output="Plink.Format" # job standard output file (%j replaced by job id)
#SBATCH --error="Plink.Format.err" # job standard error file (%j replaced by job id)

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load plink/2.0

plink2 --vcf ../Analysis_10/CA19.EqiSNPs.vcf.gz --make-bed --out CA19EqiSNPs

plink2 --vcf ../Analysis_10/CA19.ShareSNPs.vcf.gz --make-bed --out CA19.ShareSNPs

plink2 --vcf ../Analysis_10/FiltWI18.4.SharSNPs.vcf.gz --make-bed --allow-extra-chr --out FiltWI18.4.SharSNPs

plink2 --vcf ../Analysis_10/FiltWI18.5.EqiSNPs.vcf.gz --make-bed --allow-extra-chr --out FiltWI18.5.EqiSNPs

