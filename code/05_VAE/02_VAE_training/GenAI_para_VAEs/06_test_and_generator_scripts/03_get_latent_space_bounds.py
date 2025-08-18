import numpy as np
import pandas as pd
import torch
import os
import ray

for vae in [1, 2, 3, 4, 5]:
    print(f"Start VAE {vae}\n\n")
    # load model
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    encoder_path = "/mnt/Data/Dropbox/Projekte/GenAI_para_for_HLSMs/data/08_trained_vaes/vae" + str(vae) + "_encoder.pt"
    encoder = torch.load(encoder_path)
    encoder = encoder.to(device)
    encoder.eval()

    # data loader
    if device == "cuda":
        ray.init(ignore_reinit_error=True, num_cpus=20, num_gpus=2)
    else:
        ray.init(ignore_reinit_error=True, num_cpus=20)
    data_main_path = "/mnt/Data/Dropbox/Projekte/GenAI_para_for_HLSMs/data/07_VAE_data/04_vae_data/vae_" + str(vae)
    data_paths = [data_main_path + "/train_data", data_main_path + "/test_data"]
    data_files = ray.data.read_json(data_paths)

    result_path = "/mnt/Data/Dropbox/Projekte/GenAI_para_for_HLSMs/data/08_trained_vaes/"
    print(f"Test results will be saved in {result_path}")
    os.makedirs(result_path, exist_ok=True)
    max_count = 2000000
    count = 0
    maximums = np.zeros((30))
    minimums = np.zeros((30))
    batch_size = 10000

    with torch.no_grad():
        for data in data_files.iter_torch_batches(batch_size=batch_size):
            x = data["function"].to(device)
            x2 = data["quantiles"].to(device, dtype=torch.float)
            count += 1
            print("count:" + str(count))
            z, mean, log_var = encoder(x, x2)
            z_np = z.cpu().detach().numpy()
            mi = z_np.min(axis=0)
            ma = z_np.max(axis=0)
            minimums[mi < minimums] = mi[mi < minimums]
            maximums[ma > maximums] = ma[ma > maximums]
            print(f"{batch_size * count / max_count * 100:.2f}%")
            print(f"Current min: {minimums.min():.4f}, current max: {maximums.max():.4f}")
            if count * batch_size > max_count:
                break

    bounds = pd.DataFrame({"min": minimums, "max": maximums})
    bounds.to_csv(result_path + "vae" + str(vae) + "_latent_bounds.csv")
    ray.shutdown()
