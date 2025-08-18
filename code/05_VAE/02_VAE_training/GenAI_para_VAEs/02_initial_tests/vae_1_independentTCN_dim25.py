import sys
import os

vae = 1

if os.name == "nt":
    path = "/"
    data_path = "C:/Users/morit/Dropbox/Projekte/GenAI_para_for_HLSMs/data/07_VAE_data/04_vae_data/vae_" + str(vae)
else:
    path = "/gpfs/data/fs71468/GenAI_para_runs/GenAI_VAE"
    data_path = "/gpfs/data/fs71468/GenAI_para_runs/data/07_VAE_data/04_vae_data/vae_" + str(vae)

# path = "/mnt/Data/Dropbox/Projekte/GenAI_VAE"
# data_path = "/mnt/Data/Dropbox/Projekte/GenAI_para_for_HLSMs/data/07_VAE_data/04_vae_data/vae_" + str(vae)

sys.path.append(path)

from GenAI_VAE.training.trainer import train_vae
import ray

# 1. Define non optimizable parameter
non_opt_para = {
    "run_name": "vae_" + str(vae) + "_bs4_ITCN_dim25",
    "train_data": data_path + "/train_data",
    "test_data": data_path + "/test_data",
    "vocab_dir": data_path + "/vocab_obj.pth",
    "results_dir": path + "/results",
    "type": "independentTCN",
    "num_workers": 2,
    "use_gpu": True,
    "data_tune_seed": 130422,
    "val_size": 0.2,
    "epochs": 200,
    "epoch_stop": 10,
    "mse_weight": 0.0001,
    "KL_weight": 0.001,
    "latent_dim": 25,
    "max_norm": 1.0,
    "warmup_step_fraction": 0.00001,
    "cross_entropy_reduction": "mean",
    "trainer_cpus": 10,
    "scheduler_max_t": 100 * 50,  # max over number of epochs, now chosen as num_samples*50
    "scheduler_reductio_factor": 4,
    "num_samples": 100,
    "max_concurrent_trials": 2,
    "init_checkpoint": "",
    "scale": False
}

# train
config = {
    "seed": 130422,
    "learning_rate": 0.0001,
    "dropout": 0.0,
    "batch_size": 8,

    # TCN hyperparameters
    "embedding_dim": 20,
    "hidden_dim": 512,

    # encoder
    "encoder_fun_tcn_units": 512,
    "encoder_fun_hidden_dim": 512,
    "encoder_kernel_size": 8,
    "encoder_fun_tcn_layers": 8,
    # decoder
    "decoder_fun_tcn_units": 512,
    "decoder_fun_hidden_dim": 512,
    "decoder_kernel_size": 3,
    "decoder_fun_tcn_layers": 8,
    # quantiles
    "quant_tcn_units": 64,
    "quant_hidden_dim": 32,

}

# ray.init(address='auto', _redis_password=os.environ['redis_password']) multi worker version
ray.init(ignore_reinit_error=True, num_cpus=non_opt_para["trainer_cpus"] + 20, num_gpus=non_opt_para["num_workers"]) # single worker version
results = train_vae(config, non_opt_para)
