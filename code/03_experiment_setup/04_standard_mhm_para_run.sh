#!/bin/sh
#SBATCH -J standard_run
#SBATCH --nodes=1
#SBATCH --qos=zen3_0512_a100x2
#SBATCH --partition=zen3_0512_a100x2
#SBATCH --mail-user=moritz.feigl@boku.ac.at
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=04_standard_mhm_para_run.out
source ~/env/spackenv
cd /home/fs71468/mfeigl/GenAI_para_for_HLSMs/code/03_experiment_setup
Rscript 04_standard_mhm_para_run.R
