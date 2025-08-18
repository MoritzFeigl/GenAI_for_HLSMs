import os
import ray
from GenAI_VAE.training.multinode_trainer import train_vae
import os

vae = 1
path = "/gpfs/data/fs71468/GenAI_para_runs/GenAI_VAE"
data_path = "/gpfs/data/fs71468/GenAI_para_runs/data/07_VAE_data/04_vae_data/vae_" + str(vae)

ray.init(address=os.environ["ip_head"])

print("Nodes in the Ray cluster:")
print(ray.nodes())

# 1. Define non optimizable parameter
non_opt_para = {
    "run_name": "vae_" + str(vae) + "_ITCN",    # defines name of result dir: results/run_name
    "train_data": data_path + "/train_data",            # path to train data
    "test_data": data_path + "/test_data",              # path to test data
    "vocab_dir": data_path + "/vocab_obj.pth",          # path to vocabulary
    "results_dir": path + "/results",                   # path to results fir
    "type": "independentTCN",                             # Model type: TCN, standardTTCN, T2TCN
    "num_workers": 6,                                   # Number of GPU if use_gpu=True or CPU if use_gpu=False
    "use_gpu": True,                                    # Compute on GPU or CPU
    "data_tune_seed": 130422,                           # Seed for data shuffle
    "val_size": 0.5,                                    # fraction of train data used for validation
    "epochs": 100,                                      # max number of epochs
    "epoch_stop": 5,                                   # number of epochs after which training is stopped --> <3 days!
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
    "init_checkpoint": "",                              # load previously saved checkpoint (only training)
    "scale": False                                      # scale quantiles from [0, 1000] to [0, 1]
}

# train vae
config = {
    "seed": 130422,
    "learning_rate": 0.0006070540028156819,
    "dropout": 0.0,
    "batch_size": 4,
    "latent_dim": 30,

    # TCN hyperparameters
    "embedding_dim": 16,
    "hidden_dim": 256,
    # encoder
    "encoder_fun_tcn_units": 256,
    "encoder_fun_hidden_dim": 256,
    "encoder_kernel_size": 6,
    "encoder_fun_tcn_layers": 5,
    # decoder
    "decoder_fun_tcn_units": 256,
    "decoder_fun_hidden_dim": 256,
    "decoder_kernel_size": 10,
    "decoder_fun_tcn_layers": 14,
    # quantiles
    "quant_tcn_units": 81,
    "quant_hidden_dim": 51,
}
results = train_vae(config, non_opt_para)

