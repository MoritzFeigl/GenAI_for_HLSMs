#!/bin/bash
#SBATCH -J tuneTCN
#SBATCH -N 1
#SBATCH --partition=zen3_0512_a100x2
#SBATCH --qos=zen3_0512_a100x2
#SBATCH --gres=gpu:2
#SBATCH --mem=200G
#SBATCH --mail-user=moritz.feigl@boku.ac.at
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=03_tune_TCN.out

# start conda and define pythonpath for GenAI_VAE library
source /gpfs/data/fs71468/GenAI_para_runs/miniconda3/bin/activate
conda activate GenAI_VAE
export PYTHONPATH="/gpfs/data/fs71468/GenAI_para_runs/GenAI_VAE"

cd /gpfs/data/fs71468/GenAI_para_runs/GenAI_VAE/GenAI_para-ET_VAEs
/gpfs/data/fs71468/GenAI_para_runs/miniconda3/envs/GenAI_VAE/bin/python3 -u 03_tune_TCN.py > 03_tune_TCN.pyout
