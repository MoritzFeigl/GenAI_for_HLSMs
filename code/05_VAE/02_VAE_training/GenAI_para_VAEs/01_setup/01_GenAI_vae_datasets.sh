#!/bin/bash
#SBATCH -J vaedata
#SBATCH -N 1
#SBATCH --partition=zen3_0512
#SBATCH --qos=zen3_0512
#SBATCH --mail-user=moritz.feigl@boku.ac.at
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=01_GenAI_vae_datasets.out

source /gpfs/data/fs71468/GenAI_para_for_HLSMs/miniconda3/bin/activate
conda activate GenAI_VAE
export PYTHONPATH="/gpfs/data/fs71468/GenAI_para_for_HLSMs/GenAI_VAE"

cd /gpfs/data/fs71468/GenAI_para_runs/GenAI_VAE/GenAI_para_for_HLSMs/01_setup
/gpfs/data/fs71468/GenAI_para_runs/miniconda3/envs/GenAI_VAE/bin/python3 -u 01_GenAI_vae_datasets.py > 01_GenAI_vae_datasets.pyout
exit
