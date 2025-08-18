import os
import sys
path = "/gpfs/data/fs71468/GenAI_para_runs/GenAI_VAE"
sys.path.append(path)

from GenAI_VAE.utils.data import prepare_json
import glo
import random
import shutil

for vae in [1, 2, 3, 4, 5]:
    data_path = "/gpfs/data/fs71468/GenAI_para_runs/data/07_VAE_data/04_vae_data/vae_" + str(vae)

    # train data preprocess
    batch_files = glob.glob(data_path + "/vae_data_batch_*.csv")
    batch_files = [file.split(".csv")[0] for file in batch_files]
    batch_files = [file.split("/")[-1] for file in batch_files]
    prepare_json(root_dir=data_path, files=batch_files, scaling=False)


    # split into training and test
    num_test_samples = max(int(len(batch_files) * 0.2), 1)

    random.seed(130422)
    test_files = random.sample(batch_files, num_test_samples)
    os.makedirs(data_path + "/test_data", exist_ok=True)
    for file in test_files:
        old_path = data_path + "/" + file + ".json"
        new_path = data_path + "/test_data/" + file + ".json"
        shutil.move(old_path, new_path)

    os.makedirs(data_path + "/train_data", exist_ok=True)
    train_files = [file for file in batch_files if not file in test_files]
    for file in train_files:
        old_path = data_path + "/" + file + ".json"
        new_path = data_path + "/train_data/" + file + ".json"
        shutil.move(old_path, new_path)

    # move csv to folder
    os.makedirs(data_path + "/original_csv", exist_ok=True)
    for file in batch_files:
        old_path = data_path + "/" + file + ".csv"
        new_path = data_path + "/original_csv/" + file + ".csv"
        shutil.move(old_path, new_path)
