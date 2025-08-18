import sys
import os
vae = 1

path = "/gpfs/data/fs71468/GenAI_para_runs/GenAI_VAE"
data_path = "/gpfs/data/fs71468/GenAI_para_runs/GenAI_VAE/data/development_dataset_2"

path = "/mnt/Data/Dropbox/Projekte/GenAI_VAE"
data_path = "/mnt/Data/Dropbox/Projekte/GenAI_VAE/data/development_dataset_2"

sys.path.append(path + "/GenAI_VAE")

from GenAI_VAE.training.multinode_trainer import tune_vae
import ray
from ray import tune

# 1. Define non optimizable parameter
non_opt_para = {
    "run_name": "XXXX_tune_devdata2_TTCN",
    "train_data": data_path + "/train_data",
    "test_data": data_path + "/test_data",
    "vocab_dir": data_path + "/vocab_obj.pth",
    "results_dir": path + "/results",
    "type": "TCN",
    "num_workers": 2,
    "use_gpu": True,
    "data_tune_seed": 130422,
    "val_size": 0.1,
    "epochs": 200,
    "epoch_stop": 20,
    "mse_weight": 0.0001,
    "KL_weight": 0.001,
    "latent_dim": 20,
    "max_norm": 1.0,
    "warmup_step_fraction": 0.0001,
    "cross_entropy_reduction": "mean",
    "trainer_cpus": 2,  # 50,
    "scheduler_max_t": 10 * 5,  # max number of epochs, now chosen as num_samples*epochs
    "scheduler_reductio_factor": 4,
    "num_samples": 10,
    "max_concurrent_trials": 20, # for tune choose to be equal to num_workers
    "init_checkpoint": "",
    "scale": False
}

# 2. tune VAE
config = {
    "seed": tune.choice([1304]),
    "learning_rate": tune.choice([0.0001]),
    "dropout": tune.choice([0.0]),
    "batch_size": tune.choice([64]),

    # transformer hyperparameters
    "embedding_dim": tune.choice([32, 64, 128]),
    "fun_tcn_layers": tune.choice([2, 3]),
    # TCN hyperparameters
    "hidden_dim": tune.choice([256, 512, 1024]),
    "fun_tcn_units": tune.choice([256, 512, 1024]),
    "fun_hidden_dim": tune.choice([128, 256, 512]),
    "quant_tcn_units": tune.choice([64]),
    "quant_hidden_dim": tune.choice([32]),
    "kernel_size": tune.choice([3, 5, 8]),
}

# specific tune settings
tune_results = tune_vae(config, non_opt_para)
