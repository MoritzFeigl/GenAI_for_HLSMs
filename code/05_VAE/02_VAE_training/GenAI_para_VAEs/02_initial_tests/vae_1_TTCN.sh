#!/bin/bash
#SBATCH -J vae1TTCN
#SBATCH -N 1
#SBATCH --partition=zen3_0512_a100x2
#SBATCH --qos=zen3_0512_a100x2
#SBATCH --gres=gpu:2
#SBATCH --mem=300G
#SBATCH --mail-user=moritz.feigl@boku.ac.at
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=vae_1_TTCN.out

cd /gpfs/data/fs71468/GenAI_para_runs/GenAI_VAE

source /gpfs/data/fs71468/GenAI_para_runs/miniconda3/bin/activate
conda activate GenAI_VAE
export PYTHONPATH="/gpfs/data/fs71468/GenAI_para_runs/GenAI_VAE"
#redis_password=$(uuidgen)
#export redis_password

#nodes=$(scontrol show hostnames $SLURM_JOB_NODELIST) # Getting the node names
#nodes_array=( $nodes )

#node_1=${nodes_array[0]}
#ip=$node_1
#port=6379
#ip_head=$ip:$port
#export ip_head
#echo "IP Head: $ip_head"

#echo "STARTING HEAD at $node_1"
#srun --nodes=1 --ntasks=1 -w $node_1 start-head.sh $ip $redis_password &
#sleep 100

#worker_num=$(($SLURM_JOB_NUM_NODES - 1)) #number of nodes other than the head node
#for ((  i=1; i<=$worker_num; i++ ))
#do
#  node_i=${nodes_array[$i]}
#  echo "STARTING WORKER $i at $node_i"
#  srun --nodes=1 --ntasks=1 -w $node_i start-worker.sh $ip_head $redis_password &
#  sleep 10
#done

/gpfs/data/fs71468/GenAI_para_runs/miniconda3/envs/GenAI_VAE/bin/python3 -u GenAI_para-ET_VAEs/vae_1_TTCN.py > GenAI_para-ET_VAEs/vae_1_TTCN.pyout
exit
