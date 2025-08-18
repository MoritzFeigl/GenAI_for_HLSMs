import torch
import torch.nn as nn
import GenAI_VAE.models.TCN_class as TCN


class Encoder(nn.Module):

    def __init__(self, embedding_dim, fun_tcn_units, quant_tcn_units, fun_hidden_dim,
                 quant_hidden_dim, hidden_dim, latent_dim, vocab_size, max_seq_length, dropout):
        super(Encoder, self).__init__()

        # parameters
        self.quant_tcn_units = quant_tcn_units

        # activation
        self.activation = nn.SELU()

        # Function encoder
        self.embedding = nn.Embedding(vocab_size, embedding_dim)
        self.fun_tcn1 = TCN.TemporalConvNet(num_inputs=embedding_dim, num_channels=[fun_tcn_units] * 5,
                                           kernel_size=3, dropout=dropout)
        self.fun_tcn2 = TCN.TemporalConvNet(num_inputs=embedding_dim, num_channels=[fun_tcn_units] * 5,
                                           kernel_size=6, dropout=dropout)
        self.fun_tcn3 = TCN.TemporalConvNet(num_inputs=embedding_dim, num_channels=[fun_tcn_units] * 5,
                                           kernel_size=12, dropout=dropout)

        self.fun_fc1 = nn.Linear(fun_tcn_units * max_seq_length * 3, fun_hidden_dim)

        # Quantile encoder
        self.quant_tcn = TCN.TemporalConvNet(num_inputs=1, num_channels=[quant_tcn_units] * 3,
                                             kernel_size=3, dropout=dropout)
        self.quant_fc1 = nn.Linear(quant_tcn_units * 9, quant_hidden_dim)

        # Latent representation layers
        self.fc1 = nn.Linear(fun_hidden_dim + quant_hidden_dim, hidden_dim)
        self.FC_mean = nn.Linear(hidden_dim, latent_dim)
        self.FC_var = nn.Linear(hidden_dim, latent_dim)
        self.embedding_dim = embedding_dim
        self.max_seq_length = max_seq_length

    def forward(self, x, x2):
        # Function encoder
        x = self.embedding(x)
        x = torch.transpose(x, 1, 2)
        x_tcn1 = self.fun_tcn1(x)
        x_tcn2 = self.fun_tcn2(x)
        x_tcn3 = self.fun_tcn3(x)
        batch_size = x.shape[0]
        x = torch.cat((x_tcn1, x_tcn2, x_tcn3), 1)
        x = x.view(batch_size, -1)
        x = self.activation(self.fun_fc1(x))

        # Quantile encoder
        x2 = self.quant_tcn(x2.view(-1, 1, 9))
        x2 = self.activation(self.quant_fc1(x2.view(-1, self.quant_tcn_units * 9)))

        # latent representation
        x_conc = torch.cat((x, x2), 1)
        h_ = self.activation(self.fc1(x_conc))
        mean = self.FC_mean(h_)
        log_var = self.FC_var(h_)  # encoder produces mean and log of variance
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
    def __init__(self, embedding_dim, fun_tcn_units, quant_tcn_units, fun_hidden_dim,
                 quant_hidden_dim, latent_dim, vocab_size, max_seq_length, dropout):
        super(Decoder, self).__init__()
        # parameters
        self.max_seq_length = max_seq_length
        self.vocab_size = vocab_size
        self.latent_dim = latent_dim
        self.quant_tcn_units = quant_tcn_units

        # activation
        self.activation = nn.SELU()

        # Function decoder
        self.z_fun_fc = nn.Linear(1, embedding_dim)
        self.fun_tcn1 = TCN.TemporalConvNet(num_inputs=embedding_dim, num_channels=[fun_tcn_units] * 5,
                                            kernel_size=3, dropout=dropout)
        self.fun_tcn2 = TCN.TemporalConvNet(num_inputs=embedding_dim, num_channels=[fun_tcn_units] * 5,
                                            kernel_size=6, dropout=dropout)
        self.fun_tcn3 = TCN.TemporalConvNet(num_inputs=embedding_dim, num_channels=[fun_tcn_units] * 5,
                                            kernel_size=12, dropout=dropout)

        self.fun_fc1 = nn.Linear(fun_tcn_units * latent_dim * 3, fun_hidden_dim)
        self.fun_fc2 = nn.Linear(fun_hidden_dim, int(fun_hidden_dim / 2))
        self.fun_fc_out = nn.Linear(int(fun_hidden_dim / 2), max_seq_length * vocab_size)
        self.softmax = nn.LogSoftmax(dim=1)

        # Quantile decoder
        self.quant_tcn = TCN.TemporalConvNet(num_inputs=1, num_channels=[quant_tcn_units] * 3,
                                             kernel_size=3, dropout=dropout)
        self.quant_fc1 = nn.Linear(quant_tcn_units * latent_dim, quant_hidden_dim)
        self.quant_fc2 = nn.Linear(quant_hidden_dim, 9)

    def forward(self, z):
        # Function decoder
        z1 = z.view(-1, self.latent_dim, 1)
        z1 = self.activation(self.z_fun_fc(z1))
        z1 = torch.transpose(z1, 1, 2)
        z1_tcn1 = self.fun_tcn1(z1)
        z1_tcn2 = self.fun_tcn2(z1)
        z1_tcn3 = self.fun_tcn3(z1)
        batch_size = z.shape[0]
        z1 = torch.cat((z1_tcn1, z1_tcn2, z1_tcn3), 1)
        z1 = z1.view(batch_size, -1)
        z1 = self.activation(self.fun_fc1(z1))
        z1 = self.activation(self.fun_fc2(z1))
        z1 = self.fun_fc_out(z1)
        z1 = z1.view(batch_size, self.vocab_size, self.max_seq_length)
        x_hat = self.softmax(z1)

        # Quantiles decoder
        z2 = self.quant_tcn(z.view(-1, 1, self.latent_dim))
        z2 = self.activation(self.quant_fc1(z2.view(-1, self.quant_tcn_units * self.latent_dim)))
        x2_hat = self.quant_fc2(z2)
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
