#!/bin/sh
#SBATCH -J training
#SBATCH -N 1
#SBATCH --qos=zen3_0512
#SBATCH --partition=zen3_0512
#SBATCH --mail-user=moritz.feigl@boku.ac.at
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=02_training_parameter_fields.out
cd /gpfs/data/fs71468/GenAI_para_runs/code/07_analyze_results
source ~/env/spackenv
spack load /54rh3fo
spack load /da7sz6h
spack load /4b4eze3
export LD_LIBRARY_PATH=$LIBRARY_PATH
conda activate maps
python3 03_save_parameter_fields.py
