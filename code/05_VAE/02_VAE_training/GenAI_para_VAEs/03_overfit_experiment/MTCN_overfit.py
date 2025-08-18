import sys
import os

if os.name == "nt":
    path = "/"
    data_path = "/data/development_dataset_2"
else:
    path = "/gpfs/data/fs71468/GenAI_para_runs/GenAI_VAE"
    data_path = "/gpfs/data/fs71468/GenAI_para_runs/GenAI_VAE/data/development_dataset_2"
sys.path.append(path)

from GenAI_VAE.training.trainer import train_vae
import ray

# 1. Define non optimizable parameter
non_opt_para = {
    "run_name": "MTCN_overfit",    # defines name of result dir: results/run_name
    "train_data": data_path + "/train_data",            # path to train data
    "test_data": data_path + "/test_data",              # path to test data
    "vocab_dir": data_path + "/vocab_obj.pth",          # path to vocabulary
    "results_dir": path + "/results",                   # path to results fir
    "type": "MTCN",                             # Model type: TCN, standardTTCN, T2TCN
    "num_workers": 2,                                   # Number of GPU if use_gpu=True or CPU if use_gpu=False
    "use_gpu": True,                                    # Compute on GPU or CPU
    "data_tune_seed": 130422,                           # Seed for data shuffle
    "val_size": 0.1,                                    # fraction of train data used for validation
    "epochs": 200,                                      # max number of epochs
    "epoch_stop": 200,                                   # number of epochs after which training is stopped --> <3 days!
    "mse_weight": 0.00001,                               # Weight of MSE in loss function
    "KL_weight": 0.001,                                 # Weight of Kullback-Leibler divergence in loss function
    "latent_dim": 155,                                   # Size of latent space
    "max_norm": 5.0,                                    # max norm clipping for gradients -->
    "warmup_step_fraction": 0.001,                      # step fraction used as warm-up for learning rate scheduler
    "cross_entropy_reduction": "mean",                  # entropy computed as mean or sum per bawtch
    "trainer_cpus": 10,  # 50,                          # number of cpus for trainer resources besides worker
    "scheduler_max_t": 10 * 5,                          # max number of epochs, now chosen as num_samples*epochs
    "scheduler_reductio_factor": 4,                     # reduction rate of ASHA scheduler
    "num_samples": 10,                          # number of samples for hyperparameter optimization
    "max_concurrent_trials": 20,                        # for tune choose to be at least to num_workers
    "init_checkpoint": "",                              # load previously saved checkpoint (only training)
    "scale": False                                      # scale quantiles from [0, 1000] to [0, 1]
}

# train vae
config = {
    # general hyperparameter
    "seed": 130422,
    "learning_rate": 0.0001,
    "dropout": 0.0,
    "batch_size": 4,
    # MTCN hyperparameters
    "embedding_dim": 32,
    "fun_tcn_units": 512,
    "fun_hidden_dim": 2048,
    "quant_tcn_units": 128,
    "quant_hidden_dim": 64,
    "hidden_dim": 1024
}
ray.init(ignore_reinit_error=True, num_cpus=non_opt_para["trainer_cpus"] + 20, num_gpus=non_opt_para["num_workers"])
results = train_vae(config, non_opt_para)
