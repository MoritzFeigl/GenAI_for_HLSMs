#!/bin/sh
#SBATCH -J GenAI_mhm
#SBATCH -N 1
#SBATCH --qos=zen3_0512
#SBATCH --partition=zen3_0512
#SBATCH --mail-user=moritz.feigl@boku.ac.at
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=01_experiment_start.out
source ~/env/spackenv
export LD_LIBRARY_PATH=$LIBRARY_PATH
source /gpfs/data/fs71468/GenAI_para_runs/miniconda3/bin/activate
export PYTHONPATH='/gpfs/data/fs71468/GenAI_para_runs/GenAI_VAE'
conda activate GenAI_VAE
cd dummy_folder
Rscript 02_mhm_GenAI_para_runs.R
