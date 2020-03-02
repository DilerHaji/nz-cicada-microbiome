#!/bin/bash
#SBATCH --job-name=######
#SBATCH --partition=general
#SBATCH --mail-type=END
#SBATCH --mail-user=diler.haji@uconn.edu
#SBATCH -c 32

cd .

mkdir fastq
for i in */*.gz; do mv $i "fastq""/"${i%%-*}"."`basename $i`; done 

module load R

Rscript dada2.r