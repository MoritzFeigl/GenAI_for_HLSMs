import GenAI_VAE.models.TTCN as TTCN
import GenAI_VAE.models.TCN as TCN
import GenAI_VAE.models.MTCN as MTCN

def load_VAE(config):
    # Transformer for TF and TCN for quantiles
    if config["type"] == "TTCN":
        # Model initialisation
        encoder = TTCN.Encoder(config["embedding_dim"], config["nhead"], config["d_hid"],
                               config["nlayers"], config["dropout"], config["quant_tcn_units"],
                               config["latent_dim"], config["vocab_size"], config["max_seq_length"])
        decoder = TTCN.Decoder(config["embedding_dim"], config["nhead"], config["d_hid"],
                               config["dropout"], config["nlayers"], config["quant_tcn_units"],
                               config["latent_dim"], config["vocab_size"],
                               config["max_seq_length"])
        VAE = TTCN.Model(encoder=encoder, decoder=decoder)

    # TCN for TF and quantiles
    if config["type"] == "TCN":
        # Model initialisation
        if 'kernel_size' in config.keys():
                encoder = TCN.Encoder(config["embedding_dim"], config["fun_tcn_units"], config["quant_tcn_units"],
                                      config["fun_hidden_dim"], config["quant_hidden_dim"], config["hidden_dim"],
                                      config["latent_dim"], config["vocab_size"], config["max_seq_length"],
                                      config["dropout"], config["kernel_size"], config["fun_tcn_layers"])
                decoder = TCN.Decoder(config["fun_tcn_units"], config["quant_tcn_units"], config["fun_hidden_dim"],
                                      config["quant_hidden_dim"], config["latent_dim"], config["vocab_size"],
                                      config["max_seq_length"], config["dropout"],
                                      config["kernel_size"], config["fun_tcn_layers"])
        else:
            encoder = TCN.Encoder(config["embedding_dim"], config["fun_tcn_units"], config["quant_tcn_units"],
                                  config["fun_hidden_dim"], config["quant_hidden_dim"], config["hidden_dim"],
                                  config["latent_dim"], config["vocab_size"], config["max_seq_length"],
                                  config["dropout"])
            decoder = TCN.Decoder(config["fun_tcn_units"], config["quant_tcn_units"], config["fun_hidden_dim"],
                                  config["quant_hidden_dim"], config["latent_dim"], config["vocab_size"],
                                  config["max_seq_length"], config["dropout"])
        VAE = TCN.Model(encoder=encoder, decoder=decoder)

    # Transformer encoder and TCN decoder for TF and TCN for quantiles
    if config["type"] == "T2TCN":
        # Model initialisation
        encoder = TTCN.Encoder(config["embedding_dim"], config["nhead"], config["d_hid"],
                               config["nlayers"], config["dropout"], config["quant_tcn_units"],
                               config["latent_dim"], config["vocab_size"], config["max_seq_length"])

        if 'kernel_size' in config.keys():
            decoder = TCN.Decoder(config["fun_tcn_units"], config["quant_tcn_units"], config["fun_hidden_dim"],
                                  config["quant_hidden_dim"], config["latent_dim"], config["vocab_size"],
                                  config["max_seq_length"], config["dropout"], config["kernel_size"])
        else:
            decoder = TCN.Decoder(config["fun_tcn_units"], config["quant_tcn_units"], config["fun_hidden_dim"],
                                  config["quant_hidden_dim"], config["latent_dim"], config["vocab_size"],
                                  config["max_seq_length"], config["dropout"])
        VAE = TCN.Model(encoder=encoder, decoder=decoder)

    # 3 TCN with different kernel size for TF and one TCN for quantiles
    if config["type"] == "MTCN":
        # Model initialisation
        encoder = MTCN.Encoder(config["embedding_dim"], config["fun_tcn_units"], config["quant_tcn_units"],
                              config["fun_hidden_dim"], config["quant_hidden_dim"], config["hidden_dim"],
                              config["latent_dim"], config["vocab_size"], config["max_seq_length"], config["dropout"])
        decoder = MTCN.Decoder(config["embedding_dim"], config["fun_tcn_units"], config["quant_tcn_units"], config["fun_hidden_dim"],
                              config["quant_hidden_dim"], config["latent_dim"], config["vocab_size"],
                              config["max_seq_length"], config["dropout"])
        VAE = MTCN.Model(encoder=encoder, decoder=decoder)

    # TCN for TF and quantiles
    if config["type"] == "independentTCN":
        # Model initialisation
        encoder = TCN.Encoder(config["embedding_dim"], config["encoder_fun_tcn_units"], config["quant_tcn_units"],
                              config["encoder_fun_hidden_dim"], config["quant_hidden_dim"], config["hidden_dim"],
                              config["latent_dim"], config["vocab_size"], config["max_seq_length"],
                              config["dropout"], config["encoder_kernel_size"], config["encoder_fun_tcn_layers"])
        decoder = TCN.Decoder(config["decoder_fun_tcn_units"], config["quant_tcn_units"], config["decoder_fun_hidden_dim"],
                              config["quant_hidden_dim"], config["latent_dim"], config["vocab_size"],
                              config["max_seq_length"], config["dropout"],
                              config["decoder_kernel_size"], config["decoder_fun_tcn_layers"])
        VAE = TCN.Model(encoder=encoder, decoder=decoder)


    if "VAE" not in locals():
        raise ValueError('Loading VAE Model failed. Check Inputs and choose a valid VAE model type: "TCN", "TTCN", '
                         'T2TCN", "MTCN", "independentTCN".')
    return VAE