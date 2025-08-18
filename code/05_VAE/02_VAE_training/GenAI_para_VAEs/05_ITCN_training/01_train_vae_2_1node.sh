#!/bin/bash
#SBATCH -J vae2_1node
#SBATCH -N 1
#SBATCH --partition=zen3_0512_a100x2
#SBATCH --qos=zen3_0512_a100x2
#SBATCH --gres=gpu:2
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=400G
#SBATCH --mail-user=moritz.feigl@boku.ac.at
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=01_train_vae_2_1node.out

# start conda and define pythonpath for GenAI_VAE library
source /gpfs/data/fs71468/GenAI_para_runs/miniconda3/bin/activate
conda activate GenAI_VAE
export PYTHONPATH="/gpfs/data/fs71468/GenAI_para_runs/GenAI_VAE"

# Start train script
cd /gpfs/data/fs71468/GenAI_para_runs/GenAI_VAE/GenAI_para-ET_VAEs/05_ITCN_training
/gpfs/data/fs71468/GenAI_para_runs/miniconda3/envs/GenAI_VAE/bin/python3 -u 01_train_vae_2_1node.py
exit
