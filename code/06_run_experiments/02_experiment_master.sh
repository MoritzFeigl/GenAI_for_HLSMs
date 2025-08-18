#!/bin/sh
#SBATCH -J XMaster
#SBATCH -N 1
#SBATCH --mem=1G
#SBATCH -n=1
#SBATCH --cpus-per-task=1
#SBATCH --qos=zen3_0512
#SBATCH --partition=zen3_0512
#SBATCH --mail-user=moritz.feigl@boku.ac.at
#SBATCH --mail-type=FAIL
#SBATCH --output=02_experiment_master.out

cd /gpfs/data/fs71468/GenAI_para_runs/code/06_run_experiments
source ~/env/spackenv
Rscript 02_experiment_master.R
