# Load
import pandas as pd
import torch
import sys
import os

# paths, vae
path = sys.argv[1]  # current run dir
vae = sys.argv[2]
setup = sys.argv[3]
# tf = sys.argv[4]


# load model
decoder_path = "/gpfs/data/fs71468/GenAI_para_runs/data/08_trained_vaes/vae" + str(vae) + "_decoder_" + str(setup) + ".pt"
decoder = torch.load(decoder_path)
decoder.eval()

# and dir setting
storage_dir = path + "/current_vae_predictions/"
os.makedirs(storage_dir, exist_ok=True)

# while loop that produces predictions whenever new inputs are available
keep_running = True
while keep_running:
    status_file = storage_dir + "vae" + str(vae) + "_setup" + str(setup) + "_start.txt"
    file = open(status_file, mode="r")
    line = file.read(3)
    if line == "sta":
        # read inputs
        input = pd.read_csv(storage_dir + "vae" + str(vae) + "_setup" + str(setup) + "_VAEinput.csv")
        input = input.to_numpy()
        input_tensor = torch.from_numpy(input).float()
        x_hat, x2_hat = decoder(input_tensor)
        print("TF generation succeeded")
        # save softmax prediction as csv
        x_hat_np = x_hat.detach().numpy()
        x_hat_df = pd.DataFrame(x_hat_np[0, :, :])
        x_hat_df.to_csv(storage_dir + "vae" + str(vae) + "_setup" + str(setup) + "_softmax.csv", index=False)
        # save quantile prediction as csv
        x2_hat_np = x2_hat.detach().numpy()
        x2_hat_df = pd.DataFrame(x2_hat_np)
        x2_hat_df.to_csv(storage_dir + "vae" + str(vae) + "_setup" + str(setup) + "_quantiles.csv", index=False)
        print(f"Finished vae {vae} setup {setup} generation")
        # change status file to end
        with open(status_file, 'w') as file:
            file.write('end\n')
