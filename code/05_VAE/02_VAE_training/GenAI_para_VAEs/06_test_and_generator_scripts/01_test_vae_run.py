import json
from GenAI_VAE.training.test_utils import analyze_hyperopt, test_model_from_checkpoint


# Hyperparameter analysis
non_opt_para_path = ""
file1 = open(non_opt_para_path)
non_opt_para = json.load(file1)
analyze_hyperopt(non_opt_para)

# Specific model run test results
non_opt_para_path = "/gpfs/data/fs71468/GenAI_para_runs/GenAI_VAE/results/vae_1_ITCN_FINAL/non_opt_para-20230605-082659.json"
config_path = "/gpfs/data/fs71468/GenAI_para_runs/GenAI_VAE/results/vae_1_ITCN_FINAL/config20230605-082659.json"
checkpoint_path = "/gpfs/data/fs71468/GenAI_para_runs/GenAI_VAE/results/vae_1_ITCN_FINAL/TorchTrainer_f9d2c_00000_0_2023-06-05_08-26-59/checkpoint_000000"
test_data_path = "/gpfs/data/fs71468/GenAI_para_runs/data/07_VAE_data/04_vae_data/vae_1/test_data"
test_model_from_checkpoint(non_opt_para_path, config_path, checkpoint_path, test_data_path)