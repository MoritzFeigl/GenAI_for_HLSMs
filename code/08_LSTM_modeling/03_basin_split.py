import pandas as pd
import random
import os
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split


main_path = "C:/Users/morit/Dropbox/Projekte/GenAI_para_for_HLSMs/"
storage_path = main_path + "data/09_lstm_data/"
data_path = main_path + "data/09_lstm_data/"

# 1. experiment settings ------------------------------------------------------------------------------------------------
basin_attr = pd.read_csv(data_path + "static/basin_attributes.csv", sep=";")
basin_attr["basin"] = basin_attr["basin"].astype("str")
# static_attributes = [col for col in basin_attr.columns if col not in ["basin", "aspect"]]
# static_attributes = ["map", "aridity", "area", "dem", "high_prec_dur", "frac_snow"]

# sample train and validation basins
trainval_dir = main_path + f"data/06_training_validation_data/Training/"
train_val_basins = [name for name in os.listdir(trainval_dir) if os.path.isdir(os.path.join(trainval_dir, name))]
random.seed(1304)

# cluster train and validation basins
trainval_attr = basin_attr.loc[basin_attr["basin"].isin(train_val_basins)].copy()
cluster_attr = ["map", "aridity", "area", "dem", "high_prec_dur", "frac_snow"]
X = trainval_attr[cluster_attr].copy()
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
k = 5
kmeans = KMeans(n_clusters=k, random_state=1304)
trainval_attr['cluster'] = kmeans.fit_predict(X_scaled)
train_df, valid_df = train_test_split(
    trainval_attr,
    test_size=0.2,
    random_state=1304,
    stratify=trainval_attr['cluster']
)

# Extract lists of basin IDs for training and validation
train_basins = train_df['basin'].tolist()
val_basins = valid_df['basin'].tolist()

print(f"Num Train basins: {len(train_basins)}, Num Val basins: {len(val_basins)}")
# Test basins as defined in project dir
test_dir = main_path + f"data/06_training_validation_data/Validation/"
test_basins = [name for name in os.listdir(test_dir) if os.path.isdir(os.path.join(test_dir, name))]

# check if all basins have static attributes
all_basins = test_basins + train_basins + val_basins
basin_attr["basin"] = basin_attr["basin"].astype(str)
basin_check = [str(basin) not in basin_attr["basin"].tolist() for basin in all_basins]
if sum(basin_check) > 0:
    raise ValueError("Basin attributes and basin list does not fit together!")
# save
with open(main_path + "code/08_LSTM_modeling/train_basins.txt", "w") as file:
    for item in train_basins:
        file.write(f"{item}\n")
with open(main_path + "code/08_LSTM_modeling/val_basins.txt", "w") as file:
    for item in val_basins:
        file.write(f"{item}\n")
with open(main_path + "code/08_LSTM_modeling/test_basins.txt", "w") as file:
    for item in test_basins:
        file.write(f"{item}\n")