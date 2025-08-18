import torch
import torch.nn as nn
from torch.nn import functional as F


class VaeLoss(torch.nn.Module):

    def __init__(self, config):
        super(VaeLoss, self).__init__()
        self.KL_weight = config["KL_weight"]
        self.mse_weight = config["mse_weight"]
        self.cross_entropy = nn.NLLLoss(reduction=config["cross_entropy_reduction"])

    def forward(self, x_hat, x, x2_hat, x2, mean, log_var):
        """Training loss consisting of cross-entropy, MSE and KL-divergence"""
        reconstruction_loss = self.cross_entropy(x_hat, x.to(dtype=torch.long))
        kld = - 0.5 * torch.sum(1 + log_var - mean.pow(2) - log_var.exp())
        mse = F.mse_loss(x2_hat, x2)
        return reconstruction_loss + kld * self.KL_weight + mse * self.mse_weight

    def separate_losses(self, x_hat, x, x2_hat, x2, mean, log_var):
        """Function returning all losses separately"""
        reconstruction_loss = self.cross_entropy(x_hat, x.to(dtype=torch.long)).item()
        kld = - 0.5 * torch.sum(1 + log_var - mean.pow(2) - log_var.exp())
        mse = F.mse_loss(x2_hat, x2).item() * self.mse_weight
        return reconstruction_loss, kld * self.KL_weight, mse

    def mean_token_accuracy(self, x_hat, x):
        """Function returning mean token prediction accuracy"""
        index_probs = torch.exp(x_hat)
        indices = torch.argmax(index_probs, dim=1)
        all_equal_token = torch.sum(indices == x)
        numer_of_notoken = torch.min(torch.sum(indices == 0), torch.sum(x == 0))
        relevant_equal_token = all_equal_token - numer_of_notoken
        relevant_token = indices.shape[0] * indices.shape[1] - numer_of_notoken
        accuracy = relevant_equal_token / relevant_token
        return accuracy
