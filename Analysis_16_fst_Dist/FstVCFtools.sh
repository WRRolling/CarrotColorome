#!/bin/bash


#SBATCH --time=04:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=96   # 48 processor core(s) per node X 2 threads per core
#SBATCH --mem=6400M   # maximum memory per node
#SBATCH --partition=priority    # standard node(s)
#SBATCH --qos msn
#SBATCH --job-name "fst"
#SBATCH --output="fst.out" # job standard output file (%j replaced by job id)
#SBATCH --error="fst.err" # job standard error file (%j replaced by job id)


module load vcftools/0.1.16

vcftools --gzvcf ../Analysis_10/CA19.ShareSNPs.vcf.gz --weir-fst-pop WesternCA19.txt --weir-fst-pop EasternCA19.txt --out CA.fst
vcftools --gzvcf ../Analysis_10/FiltWI18.4.SharSNPs.vcf.gz --weir-fst-pop WesternWI18.txt --weir-fst-pop EasternWI18.txt --out WI.fst
