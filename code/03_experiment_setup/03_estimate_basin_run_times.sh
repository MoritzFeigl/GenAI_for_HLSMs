#!/bin/sh
#SBATCH -J run_times
#SBATCH -N 1
#SBATCH --cpus-per-task=10
#SBATCH --qos=zen3_0512
#SBATCH --partition=zen3_0512
#SBATCH --mail-user=moritz.feigl@boku.ac.at
#SBATCH --mail-type=FAIL
#SBATCH --output=03_estimate_basin_run_times.out
source ~/env/spackenv
cd /home/fs71468/mfeigl/GenAI_para_for_HLSMs/code/03_experiment_setup
Rscript 03_estimate_basin_run_times.R
