import json
import torch
from ray.air import Checkpoint
from GenAI_VAE.training.load_VAE import load_VAE
import pandas as pd
from GenAI_VAE.training.trainer import indices_to_functions
from GenAI_VAE.training.test_utils import test_model_from_checkpoint
import os

# Choos vae checkpoints with lowest val_loss
# path = "/mnt/Data/Dropbox/Projekte/GenAI_VAE/results/"
path = "/gpfs/data/fs71468/GenAI_para_runs/GenAI_VAE/results/"
vae_1 = {"vae": 1,
         "non_opt_para_path": path + "vae_1_ITCN_FINAL/non_opt_para-20230417-091558.json",
         "config_path": path + "vae_1_ITCN_FINAL/config20230417-091558.json",
         "checkpoint_path": path + "vae_1_ITCN_FINAL/TorchTrainer_b3ef5_00000_0_2023-04-17_09-15-59/checkpoint_000003",
         "data_dir": "/gpfs/data/fs71468/GenAI_para_runs/data/07_VAE_data/04_vae_data/vae_1"
         }

vae_2 = {"vae": 2,
         "non_opt_para_path": path + "vae_2_ITCN_FINAL/non_opt_para-20230417-083525.json",
         "config_path": path + "vae_2_ITCN_FINAL/config20230417-083525.json",
         "checkpoint_path": path + "vae_2_ITCN_FINAL/TorchTrainer_09a9b_00000_0_2023-04-17_08-35-26/checkpoint_000004",
         "data_dir": "/gpfs/data/fs71468/GenAI_para_runs/data/07_VAE_data/04_vae_data/vae_2"
         }

vae_3 = {"vae": 3,
         "non_opt_para_path": path + "vae_3_ITCN_1Node/non_opt_para-20230414-094043.json",
         "config_path": path + "vae_3_ITCN_1Node/config20230414-094043.json",
         "checkpoint_path": path + "vae_3_ITCN_1Node/TorchTrainer_a99b6_00000_0_2023-04-14_09-40-44/checkpoint_000009",
         "data_dir": "/gpfs/data/fs71468/GenAI_para_runs/data/07_VAE_data/04_vae_data/vae_3"
         }

vae_4 = {"vae": 4,
         "non_opt_para_path": path + "vae_4_ITCN_FINAL/non_opt_para-20230417-085830.json",
         "config_path": path + "vae_4_ITCN_FINAL/config20230417-085830.json",
         "checkpoint_path": path + "vae_4_ITCN_FINAL/TorchTrainer_42d8c_00000_0_2023-04-17_08-58-30/checkpoint_000004",
         "data_dir": "/gpfs/data/fs71468/GenAI_para_runs/data/07_VAE_data/04_vae_data/vae_4"
         }
vae_5 = {"vae": 5,
         "non_opt_para_path": path + "vae_5_ITCN_1Node/non_opt_para-20230414-093424.json",
         "config_path": path + "vae_5_ITCN_1Node/config20230414-093424.json",
         "checkpoint_path": path + "vae_5_ITCN_1Node/TorchTrainer_c7df8_00000_0_2023-04-14_09-34-25/checkpoint_000000",
         "data_dir": "/gpfs/data/fs71468/GenAI_para_runs/data/07_VAE_data/04_vae_data/vae_5"
         }


def save_vae_for_inference(vae_infos):
    # load non opt para best config
    file1 = open(vae_infos["non_opt_para_path"])
    non_opt_para = json.load(file1)

    # load config
    file2 = open(vae_infos["config_path"])
    config = json.load(file2)

    # checkpoint
    checkpoint = Checkpoint.from_directory(vae_infos["checkpoint_path"])

    # device and data directory
    device = "cpu"

    # vocab infos
    non_opt_para["vocabulary"] = torch.load(vae_infos["data_dir"] + "/vocab_obj.pth")
    non_opt_para["vocab_size"] = len(non_opt_para["vocabulary"].get_stoi())
    non_opt_para["max_seq_length"] = pd.read_csv(vae_infos["data_dir"] + "/max_seq_length.csv")["max_seq_length"][0]
    config = config | non_opt_para

    # define model & loss
    model = load_VAE(config)
    model = model.to(device)

    # load checkpoint, rename keys and load model weights
    checkpoint_dict = checkpoint.to_dict()
    weights = checkpoint_dict.get("model_weights")
    keys = list(weights.keys())
    new_keys = [k.replace("module.", "") for k in keys]

    def rename_keys(dict_, new_keys):
        """
         new_keys: type List(), must match length of dict_
        """
        # dict_ = {oldK: value}
        # d1={oldK:newK,} maps old keys to the new ones:
        d1 = dict(zip(list(dict_.keys()), new_keys))
        return {d1[oldK]: value for oldK, value in dict_.items()}

    weights = rename_keys(weights, new_keys)

    model.load_state_dict(weights)
    model.eval()

    # Specify a path
    model_storeage_dir = "/gpfs/data/fs71468/GenAI_para_runs/data/08_trained_vaes"
    os.makedirs(model_storeage_dir, exist_ok=True)

    encoder_path = model_storeage_dir + "/vae" + str(vae_infos["vae"]) + "_encoder.pt"
    decoder_path = model_storeage_dir + "/vae" + str(vae_infos["vae"]) + "_decoder.pt"

    # Save
    torch.save(model.encoder, encoder_path)
    torch.save(model.decoder, decoder_path)
    vocab = pd.DataFrame(non_opt_para["vocabulary"].get_itos())
    vocab.to_csv(model_storeage_dir + "/vae" + str(vae_infos["vae"]) + "_vocabulary.csv")

    # Load
    decoder = torch.load(decoder_path)
    decoder.eval()

    ##data = torch.normal(0, 0.5, (1, 30))
    # x_hat, x2_hat = decoder(data)
    # index_probs = torch.exp(x_hat)
    # indices = torch.argmax(index_probs, dim=1).cpu().detach().numpy()
    # predicted_tfs = indices_to_functions(indices, non_opt_para["vocabulary"].get_itos())
    # predicted_quantiles = x2_hat.cpu().detach().numpy()


for vae in [vae_1, vae_2, vae_3, vae_4, vae_5]:
    print("preparing and testing vae " + str(vae["vae"]), "\n\n")
    save_vae_for_inference(vae)
    #test_model_from_checkpoint(vae["non_opt_para_path"], vae["config_path"], vae["checkpoint_path"],
    #                           vae["data_dir"] + "/test_data")