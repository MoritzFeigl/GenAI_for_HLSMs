import sys
import os
vae = 1

path = "/gpfs/data/fs71468/GenAI_para_runs/GenAI_VAE"
data_path = "/gpfs/data/fs71468/GenAI_para_runs/GenAI_VAE/data/development_dataset_2"

# path = "/mnt/Data/Dropbox/Projekte/GenAI_VAE"
# data_path = "/mnt/Data/Dropbox/Projekte/GenAI_VAE/data/development_dataset_2"

sys.path.append(path)

from GenAI_VAE.training.trainer import tune_vae
import ray
from ray import tune

# 1. Define non optimizable parameter
non_opt_para = {
    "run_name": "tune_devdata2_MTCN",
    "train_data": data_path + "/train_data",
    "test_data": data_path + "/test_data",
    "vocab_dir": data_path + "/vocab_obj.pth",
    "results_dir": path + "/results",
    "type": "MTCN",
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
    "warmup_step_fraction": 0.001,
    "cross_entropy_reduction": "mean",
    "trainer_cpus": 10,  # 50,
    "scheduler_max_t": 10 * 5,  # max number of epochs, now chosen as num_samples*epochs
    "scheduler_reductio_factor": 4,
    "num_samples": 50,
    "max_concurrent_trials": 20, # for tune choose to be equal to num_workers
    "init_checkpoint": "", #path + "/results/vae_1_bs4_T2TCN_simple/TorchTrainer_3cbfb_00000_0_2023-02-20_11-06-26/checkpoint_000004"
    "scale": False
}

# 2. tune VAE
config = {
    "seed": tune.choice([1304]),
    "learning_rate": tune.choice([0.0001]),
    "dropout": tune.choice([0.0]),
    "batch_size": tune.choice([64]),

    # TCN hyperparameters
    "embedding_dim": tune.randint(16, 40),
    "hidden_dim": tune.choice([256, 512, 1024, 2048]),
    "fun_tcn_units": tune.choice([256, 512, 1024, 2048]),
    "fun_hidden_dim": tune.choice([256, 512, 1024, 2048]),
    "quant_tcn_units": tune.choice([64]),
    "quant_hidden_dim": tune.choice([32]),
    "kernel_size": tune.randint(2, 15)
}

# specific tune settings
#ray.init(address='auto', _redis_password=os.environ['redis_password'], include_dashboard=False)
ray.init(ignore_reinit_error=True, num_cpus=non_opt_para["trainer_cpus"] + 20, num_gpus=non_opt_para["num_workers"])
tune_results = tune_vae(config, non_opt_para)
