import pandas as pd
import torch
import os
import ray
import glob
import json
from GenAI_VAE.training.trainer import indices_to_functions

parameters = ["Ksat", "FieldCap", "ThetaS_1", "ThetaS_2", "fRoots_1", "fRoots_2", "L1_Max_Canopy_Intercept"]
ranges = {"Ksat": [5, 300],  # https://climexhandbook.w.uib.no/2019/11/05/soil-hydraulic-conductivity/ in cm/day
          "FieldCap": [0.1, 0.35], # https://en.wikipedia.org/wiki/Water_content as fraction
          "ThetaS_1": [0.20, 0.50], "ThetaS_2": [0.20, 0.50], # https://en.wikipedia.org/wiki/Water_content as fraction
          "fRoots_1": [0.1, 0.8], "fRoots_2": [0.1, 0.8], # as fraction
          "L1_Max_Canopy_Intercept": [0, 10]} # in mm

for vae in [1, 2, 3, 4, 5]:
    for parameter in parameters:
        if vae == 5 and parameter != "L1_Max_Canopy_Intercept":
            continue
        print(f"Start VAE {vae} - {parameter}\n\n")

        # generate data subset -----------------------------------------------------------------------------------------
        data_path = "/mnt/Data/Dropbox/Projekte/GenAI_para_for_HLSMs/data/07_VAE_data/04_vae_data/vae_" + str(vae)
        para_path = data_path + "/" + "parameter_specific"
        para_sub_path = para_path + "/" + parameter
        # save only entries with relevant ranges
        os.makedirs(para_path, exist_ok=True)
        os.makedirs(para_sub_path, exist_ok=True)
        for split in ["train_data", "test_data"]:
            batch_files = glob.glob(data_path + "/" + split + "/vae_data_batch_*.json")
            batch_files = [file.split(".json")[0] for file in batch_files]
            batch_files = [file.split("/")[-1] for file in batch_files]
            count = 0
            for file in batch_files:
                if count == 0:
                    data = []
                # subset train data entries
                with open(data_path + "/" + split + "/" + file + ".json") as f:
                    for line in f:
                        loaded_line = json.loads(line)
                        quantiles = loaded_line["quantiles"]
                        if quantiles[0] > ranges[parameter][0] and quantiles[8] < ranges[parameter][1]:
                            data.append(loaded_line)
                count += 1
                if count > 3 and len(data) > 100:
                    # save as json file in long_seq dir
                    with open(para_sub_path + "/" + file + ".json", 'w') as outfile:
                        for entry in data:
                            json.dump(entry, outfile)
                            outfile.write('\n')
                    count = 0

        # run encoder model --------------------------------------------------------------------------------------------
        # load model
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        encoder_path = "/mnt/Data/Dropbox/Projekte/GenAI_para_for_HLSMs/data/08_trained_vaes/vae" + str(vae) + "_encoder.pt"
        encoder = torch.load(encoder_path)
        encoder = encoder.to(device)
        encoder.eval()

        # get dataset
        if device == "cuda":
            ray.init(ignore_reinit_error=True, num_cpus=20, num_gpus=2)
        else:
            ray.init(ignore_reinit_error=True, num_cpus=20)
        data_files = ray.data.read_json(para_sub_path)

        vocabulary = pd.read_csv("/mnt/Data/Dropbox/Projekte/GenAI_para_for_HLSMs/data/08_trained_vaes/vae" + str(vae) + "_vocabulary.csv")
        vocab = vocabulary.iloc[:, 1].to_list()
        test_predictions = []
        print(f"Test results will be saved in {para_sub_path}")

        max_count = data_files.count()
        count = 0
        batch_size = 10000
        with torch.no_grad():
            predictions = []
            count = val_loss = mean_entropy = mean_token_accuracy = 0
            # get all test predictions
            for data in data_files.iter_torch_batches(batch_size=batch_size):
                x = data["function"].to(device)
                x2 = data["quantiles"].to(device, dtype=torch.float)
                count += 1
                print(f"count:{count}/{max_count/batch_size}")
                z, mean, log_var = encoder(x, x2)
                true_tfs = indices_to_functions(x.cpu(), vocab)
                true_quantiles = x2.cpu().detach().numpy()
                pred_batch = pd.concat((pd.DataFrame(true_tfs),
                                             pd.DataFrame(true_quantiles),
                                             pd.DataFrame(z.cpu())), axis=1)
                pred_batch.columns = ["function"] + \
                                          [f"Q0.{i}" for i in range(1, 10)] + [f"LS{i}" for i in range(1, 31)]
                predictions.append(pred_batch)
        # save test predictions
        pred_df = pd.concat(predictions)
        pred_df.to_csv(os.path.join(para_sub_path, "../vae"+str(vae)+"_"+parameter+"_samples.csv"), index=False)
        ray.shutdown()








