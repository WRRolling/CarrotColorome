#!/bin/bash

#SBATCH --time=12:0:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=80   # 36 processor core(s) per node X 2 threads per core
#SBATCH --mem=1494G   # maximum memory per node
#SBATCH --partition=scavenger    # standard node(s)
#SBATCH --job-name="Filter.1"
#SBATCH --mail-user=william.rolling@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output="Filter.1" # job standard output file (%j replaced by job id)
#SBATCH --error="Filter.1" # job standard error file (%j replaced by job id)

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

# load vcftools module
module load vcftools/0.1.16
module load samtools/1.17
# Subset individuals 
vcftools --gzvcf ../../Reseq.Carotenoid.Nov2020/V3.genome/Genotype.File.Processing/Raw.data/kevinG.d5.vcf.gz \
	--keep list2subset.txt \
	--recode \
	--stdout | bgzip -c > WI18.vcf.gz

# options for how to filter VCF file
vcftools --gzvcf WI18.vcf.gz  \
	--max-missing 0.3 \
	--minDP 5 \
	--remove-indels \
	--min-alleles 2 \
	--max-alleles 2 \
	--maf 0.05 \
	--max-maf 0.95 \
	--recode \
	--stdout | bgzip -c > FiltWI18.vcf.gz

