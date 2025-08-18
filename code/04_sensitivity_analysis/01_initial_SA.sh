#!/bin/sh
#SBATCH -J SA
#SBATCH -N 1
#SBATCH --qos=zen3_0512
#SBATCH --partition=zen3_0512
#SBATCH --mail-user=moritz.feigl@boku.ac.at
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=01_initial_SA.out
source ~/env/spackenv
cd /gpfs/data/fs71468/GenAI_para_runs/code/04_sensitivity_analysis
Rscript 01_initial_SA.R
