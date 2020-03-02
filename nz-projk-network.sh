#!/bin/bash
#SBATCH --job-name=micronet
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=200G
#SBATCH --mail-user=diler.haji@uconn.edu
#SBATCH -o myscript_%j.out
#SBATCH -e myscript_%j.err

cd . 

module load R/3.4.3

Rscript nzprojk-network.r
