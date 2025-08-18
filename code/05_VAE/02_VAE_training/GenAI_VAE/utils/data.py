import pandas as pd
from torch.utils.data import Dataset
import torch
import os
import json
import torchtext
from torchtext.data import functional
from collections import Counter
import numpy as np

# Currently not used, only used to prepare json data....
class TFdataset(Dataset):
    """Transfer functions dataset"""

    def __init__(self, root_dir, is_train):
        """Initializes the transfer function dataset"""
        self.root_dir = root_dir
        # get file name
        if is_train:
            self.file_name = "train_data"
        else:
            self.file_name = "test_data"
        # apply preprocessing if necessary and store as json
        if not os.path.isfile(root_dir + "/" + self.file_name + ".json"):
            self.prepare_json(root_dir=root_dir)
        else:
            # load vocabulary
            self.vocabulary = torch.load(root_dir + "/vocab_obj.pth")
        # load data
        with open(self.root_dir + "/" + self.file_name + ".json") as json_file:
            self.lines = json_file.readlines()
        # load max sequence length
        self.max_seq_length = pd.read_csv(root_dir + "/max_seq_length.csv").iloc[0, 0]

    def __len__(self):
        return len(self.lines)

    def __getitem__(self, idx):
        jdata = json.loads(self.lines[idx])
        return np.array(jdata["function"], dtype=np.int32), np.array(jdata["quantiles"], dtype=np.float32)

    def prepare_json(self, root_dir, files=["train_data", "test_data"]):
        """Prepares train_data and test_data json files with padding and min-max scaling"""
        # find max sequence length for padding
        max_seq_length = 0
        for file_name in files:
            raw_data = pd.read_csv(root_dir + "/" + file_name + ".csv")
            split_functions = list(functional.simple_space_split(list(raw_data.iloc[:, 0])))  # split on spaces
            max_seq_length = max(max_seq_length, max([len(i) for i in split_functions]))
        pd.DataFrame({'max_seq_length': [max_seq_length]}).to_csv(root_dir + "/max_seq_length.csv", index=False)

        # get min-max quantile value from training data for scaling
        raw_data = pd.read_csv(root_dir + "/train_data.csv")
        max_q = raw_data.iloc[:, 1:].max().max()
        min_q = raw_data.iloc[:, 1:].min().min()

        # Prepare train and test data as json
        for file_name in files:
            print(f"Preparing {file_name} Data Frame and saving in json")
            # load data and split functions on spaces
            raw_data = pd.read_csv(root_dir + "/" + file_name + ".csv")
            split_functions = list(functional.simple_space_split(list(raw_data.iloc[:, 0])))
            # create vocabulary
            if file_name == "train_data":
                counter = Counter([item for sublist in split_functions for item in sublist])
                self.vocabulary = torchtext.vocab.vocab(counter)
                self.vocabulary.insert_token('<pad>', 0)
                torch.save(self.vocabulary, root_dir + "/vocab_obj.pth")
            # min-max scaling of all quantile columns
            quantile_data = raw_data.iloc[:, 1:]
            quantile_data_scaled = quantile_data.apply(lambda x: self.min_max_scaling(x, min_q, max_q), axis=1)
            # store as json line data
            json_data = []
            for index, line in quantile_data_scaled.iterrows():
                n_pads = max_seq_length - len(split_functions[index])
                json_line = {}
                padded_line = split_functions[index] + ["<pad>"] * n_pads
                json_line['function'] = self.vocabulary.lookup_indices(padded_line)
                json_line['quantiles'] = list(map(float, list(line)))
                json_data.append(json_line)
            # save as json file
            with open(root_dir + "/" + file_name + ".json", 'w') as outfile:
                for entry in json_data:
                    json.dump(entry, outfile)
                    outfile.write('\n')

    @staticmethod
    def min_max_scaling(x, min_q, max_q):
        x_scaled = (x - min_q) / (max_q - min_q)
        return x_scaled


def prepare_json(root_dir, files, scaling=False):
    """Prepares train_data and test_data json files with padding and min-max scaling"""

    # helper function for fixing missing spaces
    def spacer(text):
        text = text.replace("(", ' ( ')
        return text

    # find max sequence length for padding
    max_seq_length = max_q = min_q = 0
    print("Checking max seq length and min/max quantile values\n")

    for file_name in files:
        raw_data = pd.read_csv(root_dir + "/" + file_name + ".csv")
        raw_data['function'] = raw_data['function'].apply(spacer)
        split_functions = list(functional.simple_space_split(list(raw_data.iloc[:, 0])))  # split on spaces

        # get max value were > 5% of functions exist
        lens = [len(i) for i in split_functions]
        thrsld = np.quantile(lens, 0.99)
        relevant_lengths = [le for le in lens if le < thrsld]
        max_seq_length = max(max_seq_length, max(relevant_lengths))
        print(f"Current max seq len: {max_seq_length}")
        # get min-max quantile value from training data for scaling
        raw_data = pd.read_csv(root_dir + "/" + file_name + ".csv")
        max_q = max(max_q, raw_data.iloc[:, 1:].max().max())
        min_q = min(min_q, raw_data.iloc[:, 1:].min().min())
    pd.DataFrame({'max_seq_length': [max_seq_length]}).to_csv(root_dir + "/max_seq_length.csv", index=False)
    pd.DataFrame({"max_quantile_value": [max_q],
                  "min_quantile_value": [min_q]}).to_csv(root_dir + "/quantile_value_bounds.csv", index=False)

    # Prepare train and test data as json
    for file_name in files:
        print(f"Preparing {file_name} Data Frame and saving in json")

        # load data and split functions on spaces
        raw_data = pd.read_csv(root_dir + "/" + file_name + ".csv")
        raw_data['function'] = raw_data['function'].apply(spacer)
        split_functions = list(functional.simple_space_split(list(raw_data.iloc[:, 0])))

        # create vocabulary
        if file_name == files[0]:
            counter = Counter([item for sublist in split_functions for item in sublist])
            vocabulary = torchtext.vocab.vocab(counter)
            current_tokens = vocabulary.get_itos()
            # define tokens with fixed order
            numbers = [str(t) for t in range(10)]
            parenthesis = ["(", ")"]
            general_tokens = ["<pad>", "."]
            all_tokens = general_tokens + parenthesis + numbers

            variables_pre = [tok for tok in current_tokens if tok not in all_tokens]
            # fix itos if left paranthesis mixed things up
            variables = [var for var in variables_pre if "(" not in var]
            all_tokens = all_tokens + variables
            token_dict = {all_tokens[i]: 1 for i in range(0, len(all_tokens))}
            vocabulary = torchtext.vocab.vocab(token_dict)
        else:
            counter = Counter([item for sublist in split_functions for item in sublist])
            new_vocabulary = torchtext.vocab.vocab(counter)
            current_tokens = new_vocabulary.get_itos()
            variables_pre = [tok for tok in current_tokens if tok not in all_tokens]
            new_variables = [var for var in variables_pre if "(" not in var]
            if new_variables != []:
                all_tokens = all_tokens + new_variables
                token_dict = {all_tokens[i]: 1 for i in range(0, len(all_tokens))}
                vocabulary = torchtext.vocab.vocab(token_dict)

        torch.save(vocabulary, root_dir + "/vocab_obj.pth")

        quantile_data = raw_data.iloc[:, 1:]
        # store as json line data
        json_data = []
        lengths = []
        for index, line in quantile_data.iterrows():
            n_pads = max_seq_length - len(split_functions[index])
            if n_pads < 0:
                continue
            json_line = {}
            lengths.append(len(split_functions[index]))
            padded_line = split_functions[index] + ["<pad>"] * n_pads
            json_line['function'] = vocabulary.lookup_indices(padded_line)
            json_line['quantiles'] = list(map(float, list(line)))
            json_data.append(json_line)

        # Check length distribution:
        # min-25%, 25-50% 50-75% and 75-max% quantile should be equally represented
        #quantiles = []
        #for q in [0.25, 0.5, 0.75]:
        #    quantiles.append(np.quantile(lengths, q))##

        # get nr of TFs in each quantile range
        #np.sum(lengths <= quantiles[0])
        #len(lengths) - np.sum(lengths < quantiles[0]) - np.sum(lengths > quantiles[1])
        #len(lengths) - np.sum(lengths < quantiles[1]) - np.sum(lengths > quantiles[2])
        #len(lengths) - np.sum(lengths < quantiles[2])
        # save as json file
        with open(root_dir + "/" + file_name + ".json", 'w') as outfile:
            for entry in json_data:
                json.dump(entry, outfile)
                outfile.write('\n')

        # save as json file
        with open(root_dir + "/" + file_name + ".json", 'w') as outfile:
            for entry in json_data:
                json.dump(entry, outfile)
                outfile.write('\n')


def min_max_scaling(x, min_q, max_q):
    x_scaled = (x - min_q) / (max_q - min_q)
    return x_scaled



