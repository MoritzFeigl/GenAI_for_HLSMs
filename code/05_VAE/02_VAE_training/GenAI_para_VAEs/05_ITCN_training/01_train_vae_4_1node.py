import ray
from GenAI_VAE.training.trainer import train_vae

vae = 4
path = "/gpfs/data/fs71468/GenAI_para_runs/GenAI_VAE"
data_path = "/gpfs/data/fs71468/GenAI_para_runs/data/07_VAE_data/04_vae_data/vae_" + str(vae)


# 1. Define non optimizable parameter
non_opt_para = {
    "run_name": "vae_" + str(vae) + "_ITCN_FINAL",    # defines name of result dir: results/run_name
    "train_data": data_path + "/train_data",            # path to train data
    "test_data": data_path + "/test_data",              # path to test data
    "vocab_dir": data_path + "/vocab_obj.pth",          # path to vocabulary
    "results_dir": path + "/results",                   # path to results fir
    "type": "independentTCN",                             # Model type: TCN, standardTTCN, T2TCN
    "num_workers": 2,                                   # Number of GPU if use_gpu=True or CPU if use_gpu=False
    "use_gpu": True,                                    # Compute on GPU or CPU
    "data_tune_seed": 130422,                           # Seed for data shuffle
    "val_size": 0.6,                                    # fraction of train data used for validation
    "epochs": 50,                                      # max number of epochs
    "epoch_stop": 2,                                   # number of epochs after which training is stopped --> <3 days!
    "mse_weight": 0.0001,                               # Weight of MSE in loss function
    "KL_weight": 0.001,                                 # Weight of Kullback-Leibler divergence in loss function
    "max_norm": 1.0,                                    # max norm clipping for gradients -->
    "warmup_step_fraction": 0.001,                      # step fraction used as warm-up for learning rate scheduler
    "cross_entropy_reduction": "mean",                  # entropy computed as mean or sum per batch
    "trainer_cpus": 10,  # 50,                          # number of cpus for trainer resources besides worker
    "scheduler_max_t": 10 * 5,                          # max number of epochs, now chosen as num_samples*epochs
    "scheduler_reductio_factor": 4,                     # reduction rate of ASHA scheduler
    "num_samples": 10,                          # number of samples for hyperparameter optimization
    "max_concurrent_trials": 20,                        # for tune choose to be at least to num_workers
    "init_checkpoint": "/gpfs/data/fs71468/GenAI_para_runs/GenAI_VAE/results/vae_4_ITCN_FINAL/TorchTrainer_a3de4_00000_0_2023-06-16_15-13-52/checkpoint_000000",                              # load previously saved checkpoint (only training)
    "scale": False                                      # scale quantiles from [0, 1000] to [0, 1]
}

# train vae
config = {
    "seed": 130422,
    "learning_rate": 0.0001,
    "dropout": 0.0,
    "batch_size": 2,
    "latent_dim": 30,

    # TCN hyperparameters
    "embedding_dim": 20,
    "hidden_dim": 256,
    # encoder
    "encoder_fun_tcn_units": 256,
    "encoder_fun_hidden_dim": 256,
    "encoder_kernel_size": 5,
    "encoder_fun_tcn_layers": 5,
    # decoder
    "decoder_fun_tcn_units": 1024,
    "decoder_fun_hidden_dim": 1024,
    "decoder_kernel_size": 3,
    "decoder_fun_tcn_layers": 5,
    # quantiles
    "quant_tcn_units": 80,
    "quant_hidden_dim": 64,
}

ray.init(ignore_reinit_error=True, num_cpus=non_opt_para["trainer_cpus"] + 20, num_gpus=non_opt_para["num_workers"])
print("Nodes in the Ray cluster:")
print(ray.nodes())
results = train_vae(config, non_opt_para)

