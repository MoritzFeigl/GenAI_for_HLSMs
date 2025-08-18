#!/bin/sh
#SBATCH -J tf_distribution
#SBATCH -N 1
#SBATCH --qos=zen3_0512
#SBATCH --partition=zen3_0512
#SBATCH --mail-user=moritz.feigl@boku.ac.at
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=04_estimate_distribution.out
source ~/env/spackenvcd
cd /gpfs/data/fs71468/GenAI_para_runs/code/05_VAE/01_VAE_preprocessing
Rscript 04_estimate_distributions.R
