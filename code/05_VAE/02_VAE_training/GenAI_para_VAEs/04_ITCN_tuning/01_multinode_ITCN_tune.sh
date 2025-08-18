#!/bin/bash
#SBATCH -J tuneITCN
#SBATCH -N 10
#SBATCH --partition=zen3_0512_a100x2
#SBATCH --qos=zen3_0512_a100x2
#SBATCH --gres=gpu:2
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=400G
#SBATCH --mail-user=moritz.feigl@boku.ac.at
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=01_multinode_ITCN_tune.out

# start conda and define pythonpath for GenAI_VAE library
source /gpfs/data/fs71468/GenAI_para_runs/miniconda3/bin/activate
conda activate GenAI_VAE
export PYTHONPATH="/gpfs/data/fs71468/GenAI_para_runs/GenAI_VAE"

# 1. Obtain ip adress
nodes=$(scontrol show hostnames $SLURM_JOB_NODELIST) # Getting the node names
nodes_array=( $nodes )

head_node=${nodes_array[0]}
head_node_ip=$(srun --nodes=1 --ntasks=1 -w "$head_node" hostname --ip-address)

# if we detect a space character in the head node IP, we'll
# convert it to an ipv4 address. This step is optional.
if [[ "$head_node_ip" == *" "* ]]; then
IFS=' ' read -ra ADDR <<<"$head_node_ip"
if [[ ${#ADDR[0]} -gt 16 ]]; then
  head_node_ip=${ADDR[1]}
else
  head_node_ip=${ADDR[0]}
fi
echo "IPV6 address detected. We split the IPV4 address as $head_node_ip"
fi

# 2. Starting head node
port=6379
ip_head=$head_node_ip:$port
export ip_head
echo "IP Head: $ip_head"

echo "Starting HEAD at $head_node"
srun --nodes=1 --ntasks=1 -w "$head_node" \
    ray start --head --node-ip-address="$head_node_ip" --port=$port \
    --num-cpus 10 --num-gpus 2 --block &

# 3. Start worker nodes
# optional, though may be useful in certain versions of Ray < 1.0.
sleep 10

# number of nodes other than the head node
worker_num=$((SLURM_JOB_NUM_NODES - 1))

for ((i = 1; i <= worker_num; i++)); do
    node_i=${nodes_array[$i]}
    echo "Starting WORKER $i at $node_i"
    srun --nodes=1 --ntasks=1 -w "$node_i" \
        ray start --address "$ip_head" \
        --num-cpus 10 --num-gpus 2 --block &
    sleep 5
done



# Start train script
cd /gpfs/data/fs71468/GenAI_para_runs/GenAI_VAE/GenAI_para-ET_VAEs/04_ITCN_tuning
/gpfs/data/fs71468/GenAI_para_runs/miniconda3/envs/GenAI_VAE/bin/python3 -u 01_multinode_ITCN_tune.py
exit