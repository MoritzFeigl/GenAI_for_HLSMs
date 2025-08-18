from GenAI_VAE.training.losses import VaeLoss
from GenAI_VAE.training.load_VAE import load_VAE
import os
import torch
from pathlib import Path
import numpy as np
import pandas as pd
import random
from ray.air import session, Checkpoint
from ray import train
import ray.train.torch
from ray.train.torch import TorchTrainer, TorchConfig
from ray.air.config import ScalingConfig
import ray
from ray import tune
from ray.tune.schedulers import HyperBandForBOHB
from ray.tune.search.bohb import TuneBOHB
from ray.air import RunConfig, CheckpointConfig
import json
import transformers
from ray.tune import ExperimentAnalysis
import time
from datetime import datetime

pd.set_option('display.max_colwidth', None)


class Trainer(object):
    """Trainer for GenAI_VAE VAE"""

    def __init__(self, model, optimizer, scheduler, vae_loss, config):
        self.model = model
        self.optimizer = optimizer
        self.vae_loss = vae_loss
        self.batch_size = config["batch_size"]
        self.epochs = config["epochs"]
        self.path = config["run_name"]
        self.best_valid_loss = float('inf')
        self.scheduler = scheduler
        self.max_norm = config["max_norm"]
        self.training_steps = config["training_steps"]
        self.init_checkpoint = config["init_checkpoint"]
        Path(self.path).mkdir(parents=True, exist_ok=True)

    def train_and_validate(self):

        start_epoch = 0
        checkpoint = session.get_checkpoint()
        if checkpoint:
            checkpoint_dict = checkpoint.to_dict()
            start_epoch = checkpoint_dict.get("epoch", -1) + 1
            self.model.load_state_dict(checkpoint_dict.get("model_weights"))
            self.optimizer.load_state_dict(checkpoint_dict.get("optimizer_state"))
            self.scheduler.load_state_dict(checkpoint_dict.get("scheduler_state"))

        if start_epoch == 0:
            if self.init_checkpoint != "":
                checkpoint = Checkpoint.from_directory(self.init_checkpoint)
                checkpoint_dict = checkpoint.to_dict()
                start_epoch = checkpoint_dict.get("epoch", -1) + 1
                self.model.load_state_dict(checkpoint_dict.get("model_weights"))
                self.optimizer.load_state_dict(checkpoint_dict.get("optimizer_state"))
                self.scheduler.load_state_dict(checkpoint_dict.get("scheduler_state"))

        # Fix ray bug where data is not switched to the correct device
        device = next(self.model.parameters()).device
        train_shard = session.get_dataset_shard("train")
        val_shard = session.get_dataset_shard("val")

        # Training loop
        for epoch in range(start_epoch, self.epochs):
            # train mode on
            self.model.train()
            running_loss = 0.0
            epoch_steps = 0
            mean_train_token_accuracy = 0
            for data in train_shard.iter_torch_batches(batch_size=self.batch_size):
                x = data["function"].to(device)
                x2 = data["quantiles"].to(device, dtype=torch.float)
                # zero the parameter gradients
                self.optimizer.zero_grad()
                if x.shape[0] < self.batch_size:
                    continue
                x_hat, x2_hat, mean, log_var = self.model(x, x2)
                loss = self.vae_loss(x_hat, x, x2_hat, x2, mean, log_var)
                mean_train_token_accuracy += self.vae_loss.mean_token_accuracy(x_hat, x)
                loss.backward()
                # Gradient Value Clipping
                torch.nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=self.max_norm, norm_type=2)
                self.optimizer.step()
                self.scheduler.step()
                # print statistics
                running_loss += loss.item()
                epoch_steps += 1
                if epoch_steps % 50 == 0:
                    print(f"Epoch {epoch + 1}: {epoch_steps / self.training_steps * 100:.3f}% - loss: {running_loss/epoch_steps}")

            # Validate
            val_batch_size = 128
            self.model.eval()
            count = val_loss = sum_entropy = kl_loss = mse_loss = mean_token_accuracy = 0
            with torch.no_grad():
                for data in val_shard.iter_torch_batches(batch_size=val_batch_size):
                    x = data["function"].to(device)
                    x2 = data["quantiles"].to(device, dtype=torch.float)
                    if x.shape[0] < self.batch_size:
                        continue
                    count += 1
                    x_hat, x2_hat, mean, log_var = self.model(x, x2)
                    val_batch_loss = self.vae_loss(x_hat, x, x2_hat, x2, mean, log_var)
                    entropy, kld, mse = self.vae_loss.separate_losses(x_hat, x, x2_hat, x2, mean, log_var)
                    val_loss += val_batch_loss.item()
                    sum_entropy += entropy
                    kl_loss += kld
                    mse_loss += mse
                    mean_token_accuracy += self.vae_loss.mean_token_accuracy(x_hat, x)
                    if count % 100 == 0:
                        print(
                            f"Epoch {epoch + 1}: validated samples: {count * val_batch_size} - accuracy: {mean_token_accuracy / count}")
            # save checkpoint
            state_dict = self.model.state_dict()
            optimizer_dict = self.optimizer.state_dict()
            scheduler_dict = self.scheduler.state_dict()

            checkpoint = Checkpoint.from_dict(
                dict(epoch=epoch,
                     model_weights=state_dict,
                     optimizer_state=optimizer_dict,
                     scheduler_state=scheduler_dict)
            )
            session.report({"Train_Loss": (running_loss / epoch_steps),
                            "Val_Loss": (val_loss / count),
                            "Train_Accuracy": (mean_train_token_accuracy.item() / epoch_steps),
                            "Val_Accuracy": (mean_token_accuracy.item() / count),
                            "Val_Entropy": (sum_entropy / count),
                            "Weighted_Val_MSE": (mse_loss / count),
                            "Weighted_Val_KL": (kl_loss.item() / count),
                            "Epoch": epoch + 1
                            }, checkpoint=checkpoint)


def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


def train_loop_per_worker(config, non_opt_para, train_data):

    print(config)
    # compute training steps
    samples = train_data._block_num_rows() * train_data.num_blocks()
    config["training_steps"] = np.floor(samples[0] / config["batch_size"])

    # combine trainable and non trainable parameters
    config = config | non_opt_para

    # define seed
    set_seed(config["seed"])

    # define model
    model = load_VAE(config)
    model = train.torch.prepare_model(model)

    # define loss and optimizer
    vae_loss = VaeLoss(config)
    optimizer = torch.optim.AdamW(model.parameters(), lr=config["learning_rate"])
    scheduler = transformers.get_cosine_schedule_with_warmup(optimizer,
                                                             num_warmup_steps=config["warmup_step_fraction"] *
                                                                              config["training_steps"],
                                                             num_training_steps=config["training_steps"])
    # Training
    trainer = Trainer(model, optimizer, scheduler, vae_loss, config)
    trainer.train_and_validate()


def indices_to_functions(indices, vocab_list):
    """get function string from vocabulary indices"""
    tfs = []
    for index_line in indices:
        tf = np.array(vocab_list)[list(index_line)]
        tf = [id.replace('<pad>', '') for id in tf]
        tf = ''.join(tf)
        tfs.append(tf)
    return tfs


def train_vae(config, non_opt_para):
    # ray datasets from json
    dataset = ray.data.read_json(non_opt_para["train_data"])
    train_data, val_data = dataset.train_test_split(test_size=non_opt_para["val_size"],
                                                    shuffle=True,
                                                    seed=non_opt_para["data_tune_seed"])
    train_data = train_data.repartition(non_opt_para["num_workers"])

    # store config dictionaries
    dir_name = non_opt_para["results_dir"] + "/" + non_opt_para["run_name"]
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    with open(dir_name + "/config" + datetime.now().strftime("%Y%m%d-%H%M%S") + ".json", "w") as f:
        json.dump(config, f)
    with open(dir_name + "/non_opt_para-" + datetime.now().strftime("%Y%m%d-%H%M%S") + ".json", "w") as f:
        json.dump(non_opt_para, f)

    # get vocabulary infos and sequence length
    non_opt_para["vocabulary"] = torch.load(non_opt_para["vocab_dir"])
    non_opt_para["vocab_size"] = len(non_opt_para["vocabulary"].get_stoi())
    non_opt_para["max_seq_length"] = len(train_data.take(1)[0]["function"])

    def train_vae_from_non_opt_para(config):
        train_loop_per_worker(config, non_opt_para, train_data)

    # scale configuration
    scaling_config = ScalingConfig(
        # Number of distributed workers.
        num_workers=non_opt_para["num_workers"],
        # Turn on/off GPU.
        use_gpu=non_opt_para["use_gpu"],
        # Specify resources used for trainer.
        trainer_resources={"CPU": non_opt_para["trainer_cpus"]},
        # Try to schedule workers on different nodes.
        # placement_strategy="SPREAD",
    )
    # checkpoint config
    checkpoint_config = CheckpointConfig(
        num_to_keep=2,
        checkpoint_score_attribute="Val_Loss",
        checkpoint_score_order="min"
    )
    # run cofiguration
    reporter = tune.CLIReporter(
        metric_columns=["Train_Loss", "Val_Loss", "Train_Accuracy", "Val_Accuracy", "Val_Entropy",
                        "Weighted_Val_MSE", "Weighted_Val_KL", "Epoch"]
    )

    try:
        epoch_stop = non_opt_para["epoch_stop"]
    except:
        epoch_stop = non_opt_para["epochs"]

    run_config = RunConfig(name=non_opt_para["run_name"],
                           local_dir=non_opt_para["results_dir"],
                           checkpoint_config=checkpoint_config,
                           sync_config=tune.SyncConfig(syncer=None),
                           progress_reporter=reporter,
                           stop={"training_iteration": epoch_stop})

    # trainer class with windows/linux backend for gpu support
    if os.name == "nt":
        backend = "gloo"
        trainer = TorchTrainer(
            train_vae_from_non_opt_para,
            train_loop_config=config,
            scaling_config=scaling_config,
            run_config=run_config,
            datasets={"train": train_data, "val": val_data},
            torch_config=TorchConfig(backend=backend)
        )
    else:
        trainer = TorchTrainer(
            train_vae_from_non_opt_para,
            train_loop_config=config,
            scaling_config=scaling_config,
            run_config=run_config,
            datasets={"train": train_data, "val": val_data}
        )

    results = trainer.fit()
    return results


def tune_vae(config, non_opt_para):
    # ray datasets from json
    dataset = ray.data.read_json(non_opt_para["train_data"])

    # scale
    if 'scale' in non_opt_para.keys():
        if non_opt_para["scale"]:
            def transform_batch(df: list) -> list:
                df["quantiles"] = [i / 1000 for i in df["quantiles"]]
                return df

            dataset = dataset.map_batches(transform_batch)

    train_data, val_data = dataset.train_test_split(test_size=non_opt_para["val_size"],
                                                    shuffle=True,
                                                    seed=non_opt_para["data_tune_seed"])
    train_data = train_data.repartition(non_opt_para["num_workers"])
    test_data = ray.data.read_json(non_opt_para["test_data"])

    # store config dictionaries
    dir_name = non_opt_para["results_dir"] + "/" + non_opt_para["run_name"]
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    with open(dir_name + "/non_opt_para-" + datetime.now().strftime("%Y%m%d-%H%M%S") + ".json", "w") as f:
        json.dump(non_opt_para, f)

    # get vocabulary infos and sequence length
    non_opt_para["vocabulary"] = torch.load(non_opt_para["vocab_dir"])
    non_opt_para["vocab_size"] = len(non_opt_para["vocabulary"].get_stoi())
    non_opt_para["max_seq_length"] = len(train_data.take(1)[0]["function"])

    # train function definition
    def train_vae_from_non_opt_para(config):
        train_loop_per_worker(config, non_opt_para, train_data)

    # ray tune and trainer configuration
    scaling_config = ScalingConfig(
        num_workers=non_opt_para["num_workers"],
        use_gpu=non_opt_para["use_gpu"],
        resources_per_worker={"GPU": 1, "CPU": 1},
        trainer_resources={"CPU": non_opt_para["trainer_cpus"]},
        placement_strategy="SPREAD",  # Try to schedule workers on different nodes.
        _max_cpu_fraction_per_node=0.8
    )
    checkpoint_config = CheckpointConfig(num_to_keep=1,
                                         checkpoint_score_attribute="Val_Loss",
                                         checkpoint_score_order="min")
    reporter = tune.CLIReporter(
        metric_columns=["Train_Loss", "Val_Loss", "Train_Accuracy", "Val_Accuracy", "Val_Entropy",
                        "Weighted_Val_MSE", "Weighted_Val_KL", "Epoch"]
    )

    run_config = RunConfig(name=non_opt_para["run_name"],
                           local_dir=non_opt_para["results_dir"],
                           checkpoint_config=checkpoint_config,
                           # sync_config=tune.SyncConfig(syncer=None),
                           progress_reporter=reporter)

    bohb_hyperband = HyperBandForBOHB(
        time_attr="training_iteration",
        max_t=non_opt_para["scheduler_max_t"],
        reduction_factor=non_opt_para["scheduler_reductio_factor"],
        stop_last_trials=False,
    )
    bohb_search = TuneBOHB(seed=non_opt_para["data_tune_seed"])
    bohb_search = tune.search.ConcurrencyLimiter(bohb_search, max_concurrent=non_opt_para["max_concurrent_trials"])

    tune_config = tune.TuneConfig(
        metric="Val_Loss",
        mode="min",
        scheduler=bohb_hyperband,
        search_alg=bohb_search,
        num_samples=non_opt_para["num_samples"],
    )
    # trainer class
    trainer = TorchTrainer(
        train_vae_from_non_opt_para,
        scaling_config=scaling_config,
        run_config=run_config,
        datasets={"train": train_data, "val": val_data}
    )

    # function for defining trial name
    def trial_str_creator(trial):
        trialname = non_opt_para["run_name"] + "_" + non_opt_para["type"]
        return trialname + "_" + time.strftime("%Y-%m-%d_%H:%M-%S_") + trial.trial_id

    # if os.path.exists(non_opt_para["results_dir"] + "/" + non_opt_para["run_name"]):
    #     tuner = tune.Tuner.restore(
    #         path=non_opt_para["results_dir"] + "/" + non_opt_para["run_name"]
    #     )
    # else:

    tuner = tune.Tuner(
        trainer,
        param_space={"train_loop_config": config},
        tune_config=tune_config,
        _tuner_kwargs={"trial_name_creator": trial_str_creator}
    )

    set_seed(non_opt_para["data_tune_seed"])
    results = tuner.fit()

    # save best results
    result_path = non_opt_para["results_dir"] + "/" + non_opt_para["run_name"]
    best_result = results.get_best_result()
    save_best = {"config": best_result.config,
                 "metrics": best_result.metrics,
                 "best_checkpoint": str(best_result.log_dir)}

    with open(result_path + "/best_result.json", "w") as outfile:
        json.dump(save_best, outfile)
    # save data frame results
    analysis = ExperimentAnalysis(result_path)
    analysis.dataframe().to_excel(result_path + "/hyperopt_results.xlsx")

    # get best parameter set
    best_result = results.get_best_result("Val_Loss", "min")
    print("Best trial config: {}".format(best_result.config["train_loop_config"]))
    print(f"Best trial final validation accuracy: {best_result.metrics['Val_Accuracy']:.3f}")
    print(f"Best trial final validation loss: {best_result.metrics['Val_Loss']:.3f}")
    print(f"Best trial final validation entropy: {best_result.metrics['Val_Entropy']:.3f}")
    print(f"Best trial final validation MSE: {best_result.metrics['Weighted_Val_MSE']:.3f}")
    print(f"Best trial final validation KL: {best_result.metrics['Weighted_Val_KL']:.3f}")

    return results

