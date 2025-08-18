import GenAI_VAE.models.TCN_class as TCN
import GenAI_VAE.models.transformer as transformer
from torch import nn
import torch


class Encoder(nn.Module):

    def __init__(self, embedding_dim, nhead, d_hid, nlayers, dropout, quant_tcn_units,
                  latent_dim, vocab_size, max_seq_length):
        super(Encoder, self).__init__()

        # parameters
        # quant_hidden_dim, hidden_dim,fun_hidden_dim
        self.quant_tcn_units = [quant_tcn_units] * 3
        self.activation = nn.GELU()

        # Function encoder
        self.fun_encoder = transformer.TransformerEncoderModel(vocab_size, embedding_dim, nhead, d_hid, nlayers,
                                                               dropout)
        # self.fun_fc1 = nn.Linear(embedding_dim * max_seq_length, fun_hidden_dim)

        # Quantile encoder
        self.quant_tcn = TCN.TemporalConvNet(num_inputs=1, num_channels=self.quant_tcn_units,
                                             kernel_size=3, dropout=dropout, activation=self.activation)
        # self.quant_fc1 = nn.Linear(self.quant_tcn_units[0] * 9, quant_hidden_dim)

        # Latent representation layers
        #self.fc1 = nn.Linear(fun_hidden_dim + quant_hidden_dim, hidden_dim)
        self.FC_mean = nn.Linear(max_seq_length * embedding_dim + 9 * quant_tcn_units, latent_dim)
        self.FC_var = nn.Linear(max_seq_length * embedding_dim + 9 * quant_tcn_units, latent_dim)
        self.embedding_dim = embedding_dim
        self.max_seq_length = max_seq_length

    def forward(self, x, x2):
        # Function encoder
        x = self.fun_encoder(x)
        batch_size = x.shape[0]
        x = x.view(batch_size, -1)
        # x = self.activation(self.fun_fc1(x))

        # Quantile encoder
        #x2 = x2.to(dtype=torch.float32)
        x2 = self.quant_tcn(x2.view(-1, 1, 9))
        # x2 = self.activation(self.quant_fc1(x2.view(-1, self.quant_tcn_units[0] * 9)))
        x2 = x2.view(-1, self.quant_tcn_units[0] * 9)

        # latent representation
        x_conc = torch.cat((x, x2), 1)
        # x_conc = self.activation(self.fc1(x_conc))
        mean = self.FC_mean(x_conc)
        log_var = self.FC_var(x_conc)  # encoder produces mean and log of variance
        #             (i.e., parateters of simple tractable normal distribution "q")
        std = torch.exp(0.5 * log_var)  # takes exponential function
        z = self.reparameterization(mean, std)

        return z, mean, log_var

    @staticmethod
    def reparameterization(mean, std):
        epsilon = torch.rand_like(std)  # sampling epsilon
        z = mean + (std * epsilon)  # reparameterization trick
        return z


class Decoder(nn.Module):
    def __init__(self, embedding_dim, nhead, d_hid, dropout, nlayers, quant_tcn_units,
                 latent_dim, vocab_size, max_seq_length):
        super(Decoder, self).__init__()
        # parameters
        self.max_seq_length = max_seq_length
        self.vocab_size = vocab_size
        self.latent_dim = latent_dim
        self.quant_tcn_units = [quant_tcn_units] * 3

        # activation
        self.activation = nn.GELU()

        # Function decoder
        decoder_layers = nn.TransformerEncoderLayer(embedding_dim, nhead, d_hid, dropout, activation=self.activation)
        self.fun_decoder = nn.TransformerEncoder(decoder_layers, nlayers)
        self.z_fun_fc = nn.Linear(latent_dim, embedding_dim)
        self.fun_fc_out = nn.Linear(embedding_dim, max_seq_length * vocab_size)
        self.softmax = nn.LogSoftmax(dim=1)

        # Quantile decoder
        self.quant_tcn = TCN.TemporalConvNet(num_inputs=1, num_channels=self.quant_tcn_units,
                                             kernel_size=3, dropout=dropout,
                                             activation=self.activation)
        self.quant_fc = nn.Linear(quant_tcn_units * latent_dim, 9)

    def forward(self, z):
        # Function decoder
        z1 = self.activation(self.z_fun_fc(z))
        z1 = self.fun_decoder(z1)
        z1 = self.fun_fc_out(z1)
        z1 = z1.view(-1, self.vocab_size, self.max_seq_length)
        x_hat = self.softmax(z1)

        # Quantiles decoder
        z2 = self.quant_tcn(z.view(-1, 1, self.latent_dim))
        z2 = z2.view(-1, self.quant_tcn_units[0] * self.latent_dim)
        x2_hat = self.quant_fc(z2)
        return x_hat, x2_hat


class Model(nn.Module):
    def __init__(self, encoder, decoder):
        super(Model, self).__init__()
        self.encoder = encoder
        self.decoder = decoder

    def forward(self, x, x2):
        z, mean, log_var = self.encoder(x, x2)
        x_hat, x2_hat = self.decoder(z)

        return x_hat, x2_hat, mean, log_var