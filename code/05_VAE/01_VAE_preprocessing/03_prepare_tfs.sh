#!/bin/sh
#SBATCH -J prepare_tfs
#SBATCH -N 1
#SBATCH --qos=zen3_0512
#SBATCH --partition=zen3_0512
#SBATCH --mail-user=moritz.feigl@boku.ac.at
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=03_prepare_tfs.out
source ~/env/spackenvcd
cd /gpfs/data/fs71468/GenAI_para_runs/code/05_VAE/01_VAE_preprocessing
Rscript 03_prepare_tfs.R
