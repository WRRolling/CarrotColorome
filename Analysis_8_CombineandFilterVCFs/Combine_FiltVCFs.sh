#!/bin/bash

#SBATCH --time=12:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=80   # 36 processor core(s) per node X 2 threads per core
#SBATCH --mem=1494G   # maximum memory per node
#SBATCH --partition=priority-mem    # large-memory node(s)
#SBATCH --qos msn-mem
#SBATCH --job-name="Combine and Filter"
#SBATCH --mail-user=william.rolling@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output="Comb_Filt.Out" # job standard output file (%j replaced by job id)
#SBATCH --error="Comb_Filt.Error" # job standard error file (%j replaced by job id)

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

# load vcftools module
module load vcftools/0.1.16
module load samtools/1.17
module load bcftools


#vcftools --gzvcf ../../Reseq.Carotenoid.Nov2020/V3.genome/Genotype.File.Processing/Raw.data/kevinG.d5.vcf.gz \
        --keep list2subset.txt \
        --recode \
        --stdout | bgzip -c > WI18.vcf.gz

#  Initial Filter of Resequencing data
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

#Change name of chromosomes to merge with GBS data!
bcftools annotate \
	--rename-chrs chr_name_conv.txt \
	 -o FiltWI18Chr.vcf.gz \
	-Oz FiltWI18.vcf.gz

 #Index Vcfs
 tabix -p vcf  FiltWI18Chr.vcf.gz
 tabix -p vcf scott2.filtered.6.vcf.gz

#Combine vcfs
bcftools merge \
	-m all \
	-O z \
	-o GenotypeFile.vcf.gz \
	FiltWI18Chr.vcf.gz scott2.filtered.6.vcf.gz


vcftools --gzvcf GenotypeFile.vcf.gz \
        --max-missing 0.65 \
        --minDP 5 \
        --remove-indels \
        --min-alleles 2 \
        --max-alleles 2 \
        --maf 0.05 \
        --max-maf 0.95 \
        --recode \
        --stdout | bgzip -c > GenotypeFile.F1.vcf.gz


perl ../../bb.vcf --infile="GenotypeFile.F1.vcf.gz" --outfile="GenotypeFile.F2.vcf.gz" --task=het --minratio=0.3 --maxratio=0.7


vcftools --gzvcf GenotypeFile.F2.vcf.gz \
        --max-missing 0.3 \
        --minDP 5 \
        --remove-indels \
        --min-alleles 2 \
        --max-alleles 2 \
        --maf 0.05 \
        --max-maf 0.95 \
        --recode \
        --stdout | bgzip -c > GenotypeFile.F3.vcf.gz 
