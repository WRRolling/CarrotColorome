#!/bin/bash


#SBATCH --time=04:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=40   # 48 processor core(s) per node X 2 threads per core
#SBATCH --mem=6400M   # maximum memory per node
#SBATCH --qos msn
#SBATCH --partition=priority    # standard node(s)
#SBATCH --job-name="fastStructure"
#SBATCH --output="fastStructure.out" # job standard output file (%j replaced by job id)
#SBATCH --error="fastStructure.error" # job standard error file (%j replaced by job id)

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load faststructure/1.0

structure.py -K 1 --input=plinkGenos/GenotypeFile --output=fastStructure/CA.1
structure.py -K 2 --input=plinkGenos/GenotypeFile --output=fastStructure/CA.2
structure.py -K 3 --input=plinkGenos/GenotypeFile --output=fastStructure/CA.3
structure.py -K 4 --input=plinkGenos/GenotypeFile --output=fastStructure/CA.4
structure.py -K 5 --input=plinkGenos/GenotypeFile --output=fastStructure/CA.5
structure.py -K 6 --input=plinkGenos/GenotypeFile --output=fastStructure/CA.6
structure.py -K 7 --input=plinkGenos/GenotypeFile --output=fastStructure/CA.7
structure.py -K 8 --input=plinkGenos/GenotypeFile --output=fastStructure/CA.8
structure.py -K 9 --input=plinkGenos/GenotypeFile --output=fastStructure/CA.9
structure.py -K 10 --input=plinkGenos/GenotypeFile --output=fastStructure/CA.10
