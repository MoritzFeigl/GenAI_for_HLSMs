from torch.utils.data.distributed import DistributedSampler
from GenAI_VAE.utils.data import TFdataset
from torch.utils.data import DataLoader, random_split
import torch


# Load datasets
def create_data_loader(data_dir, batch_size):
    full_training_data = TFdataset(root_dir=data_dir, is_train=True)
    print(int(full_training_data.__len__()))
    lengths = [int(full_training_data.__len__() * 0.8), int(full_training_data.__len__() * 0.2)]
    train_data, val_data = random_split(full_training_data, lengths, generator=torch.Generator().manual_seed(250690))
    test_data = TFdataset(root_dir=data_dir, is_train=False)

    # load iterators
    train_sampler = DistributedSampler(dataset=train_data, shuffle=True)
    train_iterator = DataLoader(train_data, batch_size=batch_size, shuffle=True,
                                sampler=train_sampler, num_workers=10, pin_memory=True)

    val_sampler = DistributedSampler(dataset=val_data, shuffle=True)
    val_iterator = DataLoader(val_data, batch_size=batch_size, shuffle=True,
                              sampler=val_sampler, num_workers=10, pin_memory=True)

    test_sampler = DistributedSampler(dataset=test_data, shuffle=True)
    test_iterator = DataLoader(test_data, batch_size=batch_size, shuffle=False,
                               sampler=test_sampler, num_workers=10, pin_memory=True)

    # define sequence and vocab properties
    max_seq_length = test_data.max_seq_length
    vocab_size = len(full_training_data.vocabulary.get_stoi())
    return train_iterator, train_sampler, val_iterator, val_sampler, test_iterator, test_sampler, max_seq_length, vocab_size
