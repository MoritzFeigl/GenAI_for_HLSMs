#!/bin/sh
#SBATCH -J val1km
#SBATCH -N 1
#SBATCH --qos=zen3_0512
#SBATCH --partition=zen3_0512
#SBATCH --mail-user=moritz.feigl@boku.ac.at
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=01_analyze_validation_1km.out
cd /gpfs/data/fs71468/GenAI_para_runs/code/07_analyze_results
source ~/env/spackenv
Rscript 01_analyze_validation_1km.R
