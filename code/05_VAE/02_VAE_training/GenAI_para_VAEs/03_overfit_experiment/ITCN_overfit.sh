#!/bin/bash
#SBATCH -J oITCN
#SBATCH -N 1
#SBATCH --partition=zen3_0512_a100x2
#SBATCH --qos=zen3_0512_a100x2
#SBATCH --gres=gpu:2
#SBATCH --mem=200G
#SBATCH --mail-user=moritz.feigl@boku.ac.at
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=ITCN_overfit.out



source /gpfs/data/fs71468/GenAI_para_runs/miniconda3/bin/activate
conda activate GenAI_VAE
export PYTHONPATH="/gpfs/data/fs71468/GenAI_para_runs/GenAI_VAE"
cd /gpfs/data/fs71468/GenAI_para_runs/GenAI_VAE/GenAI_para-ET_VAEs/experiments/overfit_experiment
/gpfs/data/fs71468/GenAI_para_runs/miniconda3/envs/GenAI_VAE/bin/python3 -u ITCN_overfit.py > ITCN_overfit.pyout
exit
