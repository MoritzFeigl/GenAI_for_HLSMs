import matplotlib.pyplot as plt
import os
import torch
import seaborn as sns
from ray import tune
import ray
import plotly.express as px
import pandas as pd
import json
from GenAI_VAE.training.load_VAE import load_VAE
from GenAI_VAE.training.trainer import indices_to_functions
from GenAI_VAE.training.losses import VaeLoss
from ray.air import Checkpoint


def analyze_hyperopt(non_opt_para):
    # Analyse HyperOpt results
    hyperopt_path = non_opt_para["results_dir"] + "/" + non_opt_para["run_name"]
    # create hyperopt analysis dir
    analysis_dir = hyperopt_path + "/hyperopt_analysis"
    os.makedirs(analysis_dir, exist_ok=True)

    print(f"All analysis result plots are stored in {analysis_dir}")
    ray.init()
    results = tune.ExperimentAnalysis(experiment_checkpoint_path=hyperopt_path)
    results_df = results.results_df
    # fix hyperparameter names
    df_names = results_df.columns
    results_df.columns = [txt.split("/")[-1] for txt in df_names]
    hyper_para_names = results_df.loc[:, "learning_rate":].columns
    # plot hyperopt results as parallel coordinates
    fig = px.parallel_coordinates(results_df,
                                  color="Val_Accuracy",
                                  dimensions=hyper_para_names,
                                  color_continuous_scale=px.colors.diverging.Tealrose_r,
                                  # color_continuous_midpoint=0.5,
                                  width=1000, height=500)
    fig.write_image(analysis_dir + "/" + "/hyperopt_pc_plot.png", scale=2)

    # scatter plots
    for para in hyper_para_names:
        fig2 = px.scatter(results_df, x=para, y="Val_Accuracy", log_x=para == "learning_rate",
                          labels={para: para,
                                  "Val_Accuracy": "Validation accuracy",
                                  "Train_Loss": "Loss"},
                          color="Train_Loss", color_continuous_scale=px.colors.diverging.Tealrose
                          )
        fig2.write_image(analysis_dir + "/" + para + "_opt_results.png", scale=2)
    ray.shutdown()


def test_best_model(non_opt_para):

    # load hyperopt best config
    hyperopt_path = non_opt_para["results_dir"] + "/" + non_opt_para["run_name"]
    f = open(hyperopt_path + '/best_result.json')
    best_result = json.load(f)
    config = best_result["config"]["train_loop_config"]

    # best checkpoint
    ray.init()
    results = tune.ExperimentAnalysis(experiment_checkpoint_path=hyperopt_path)
    best_trial = results.get_best_trial(metric="Val_Loss", mode="min")
    checkpoint = results.get_best_checkpoint(best_trial, metric="Val_Loss", mode="min")
    #config = results.get_all_configs(metric="Val_Loss", mode="min")

    # device and data directory
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    # data loader
    test_data = ray.data.read_json(non_opt_para["test_data"])

    # vocab infos
    non_opt_para["vocabulary"] = torch.load(non_opt_para["vocab_dir"])
    non_opt_para["vocab_size"] = len(non_opt_para["vocabulary"].get_stoi())
    non_opt_para["max_seq_length"] = len(test_data.take(1)[0]["function"])
    config = config | non_opt_para

    # define model & loss
    model = load_VAE(config)
    model = model.to(device)
    vae_loss = VaeLoss(config)

    # load checkpoint
    checkpoint_dict = checkpoint.to_dict()
    model.load_state_dict(checkpoint_dict.get("model_weights"))
    test_result_path = hyperopt_path + "/best_model_test_results"
    print(f"Test results will be saved in {test_result_path}")
    os.makedirs(test_result_path, exist_ok=True)

    model.eval()
    with torch.no_grad():
        test_predictions = []
        count = val_loss = mean_entropy = mean_token_accuracy = 0
        # get all test predictions
        for data in test_data.iter_torch_batches(batch_size=config["batch_size"]):
            x = data["function"].to(device)
            x2 = data["quantiles"].to(device, dtype=torch.float)
            if x.shape[0] < config["batch_size"]:
                continue
            count += 1
            x_hat, x2_hat, mean, log_var = model(x, x2)
            index_probs = torch.exp(x_hat)
            indices = torch.argmax(index_probs, dim=1).cpu().detach().numpy()
            true_tfs = indices_to_functions(x.cpu(), non_opt_para["vocabulary"].get_itos())
            predicted_tfs = indices_to_functions(indices, non_opt_para["vocabulary"].get_itos())
            true_quantiles = x2.cpu().detach().numpy()
            predicted_quantiles = x2_hat.cpu().detach().numpy()
            test_pred_batch = pd.concat((pd.DataFrame(true_tfs),
                                         pd.DataFrame(predicted_tfs),
                                         pd.DataFrame(predicted_quantiles),
                                         pd.DataFrame(true_quantiles)), axis=1)
            test_pred_batch.columns = ["function", "predicted_function"] + \
                                      [f"pred_Q0.{i}" for i in range(1, 10)] + [f"Q0.{i}" for i in range(1, 10)]
            test_predictions.append(test_pred_batch)
            val_loss += vae_loss(x_hat, x, x2_hat, x2, mean, log_var)
            mean_token_accuracy += vae_loss.mean_token_accuracy(x_hat, x)
        # save test predictions
        test_pred_df = pd.concat(test_predictions)
        test_pred_df.to_csv(os.path.join(test_result_path, "test_predictions.csv"), index=False)
        # store test losses
        test_losses = pd.DataFrame({"test_loss": [val_loss.cpu().detach().numpy() / count],
                                    "mean_token_accuracy": [mean_token_accuracy.cpu().detach().numpy() / count]})
        test_losses.to_csv(os.path.join(test_result_path, "test_losses.csv"), index=False)

        print(f"Testing Loss: {val_loss / count:.4f}")
        print(f"Average token reconstruction accuracy: {mean_token_accuracy / count:.4f}")

        # Plot predicted vs. true quantiles of 10 random functions
        number_of_samples = 10
        test_sample = test_pred_df.sample(number_of_samples, random_state=42)  # .transpose()
        melted_samples = []
        for col in test_sample.columns:
            # only melt quantile columns
            if test_sample[col].dtype != 'float32':
                continue
            col_sub = pd.DataFrame({"function": range(1, (number_of_samples + 1)),
                                    "quantile": float(col.split("Q")[1]),
                                    "parameter": test_sample[col]})
            # define split column
            if "pred" in col:
                col_sub["split"] = "prediction"
            else:
                col_sub["split"] = "true"
            melted_samples.append(col_sub)

        samples_plt_df = pd.concat(melted_samples).reset_index()
        palette = sns.color_palette("tab10", number_of_samples)
        sns.set_style("whitegrid")
        sns.set(rc={'figure.figsize': (14, 12)})
        sns.lineplot(data=samples_plt_df,
                     y="quantile",
                     x="parameter",
                     hue="function",
                     style="split",
                     style_order=["true", "prediction"],
                     palette=palette)
        plt.savefig(os.path.join(test_result_path, "quantile_sample_plot.png"))
        test_sample.to_csv(os.path.join(test_result_path, "quantile_sample_plot.csv"), index=False)
    ray.shutdown()

def test_model_from_checkpoint(non_opt_para_path, config_path, checkpoint_path, test_data_path):

    # load non opt para best config
    file1 = open(non_opt_para_path)
    non_opt_para = json.load(file1)

    # load config
    file2 = open(config_path)
    config = json.load(file2)

    # checkpoint
    checkpoint = Checkpoint.from_directory(checkpoint_path)

    # device and data directory
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    # data loader
    if device == "cuda":
        ray.init(ignore_reinit_error=True, num_cpus=20, num_gpus=non_opt_para["num_workers"])
    else:
        ray.init(ignore_reinit_error=True, num_cpus=20)
    test_data = ray.data.read_json(test_data_path)

    # vocab infos
    non_opt_para["vocabulary"] = torch.load(non_opt_para["vocab_dir"])
    non_opt_para["vocab_size"] = len(non_opt_para["vocabulary"].get_stoi())
    non_opt_para["max_seq_length"] = len(test_data.take(1)[0]["function"])
    config = config | non_opt_para

    # define model & loss
    model = load_VAE(config)
    model = torch.nn.DataParallel(model)
    model = model.to(device)
    vae_loss = VaeLoss(config)

    # load checkpoint
    checkpoint_dict = checkpoint.to_dict()
    model.load_state_dict(checkpoint_dict.get("model_weights"))
    test_result_path = non_opt_para["results_dir"] + "/" + non_opt_para["run_name"] + "/test_results"
    print(f"Test results will be saved in {test_result_path}")
    os.makedirs(test_result_path, exist_ok=True)

    model.eval()
    with torch.no_grad():
        test_predictions = []
        count = val_loss = mean_entropy = mean_token_accuracy = 0
        # get all test predictions
        for data in test_data.iter_torch_batches(batch_size=1024):#config["batch_size"]):
            x = data["function"].to(device)
            x2 = data["quantiles"].to(device, dtype=torch.float)
            if x.shape[0] < config["batch_size"]:
                continue
            count += 1
            print("count:" + str(count))
            x_hat, x2_hat, mean, log_var = model(x, x2)
            index_probs = torch.exp(x_hat)
            indices = torch.argmax(index_probs, dim=1).cpu().detach().numpy()
            true_tfs = indices_to_functions(x.cpu(), non_opt_para["vocabulary"].get_itos())
            predicted_tfs = indices_to_functions(indices, non_opt_para["vocabulary"].get_itos())
            true_quantiles = x2.cpu().detach().numpy()
            predicted_quantiles = x2_hat.cpu().detach().numpy()
            test_pred_batch = pd.concat((pd.DataFrame(true_tfs),
                                         pd.DataFrame(predicted_tfs),
                                         pd.DataFrame(predicted_quantiles),
                                         pd.DataFrame(true_quantiles)), axis=1)
            test_pred_batch.columns = ["function", "predicted_function"] + \
                                      [f"pred_Q0.{i}" for i in range(1, 10)] + [f"Q0.{i}" for i in range(1, 10)]
            test_predictions.append(test_pred_batch)
            val_loss += vae_loss(x_hat, x, x2_hat, x2, mean, log_var)
            mean_token_accuracy += vae_loss.mean_token_accuracy(x_hat, x)
            print(f"Mean acc: {mean_token_accuracy/count}")
        # save test predictions
        test_pred_df = pd.concat(test_predictions)
        test_pred_df.to_csv(os.path.join(test_result_path, "test_predictions.csv"), index=False)
        # store test losses
        test_losses = pd.DataFrame({"test_loss": [val_loss.cpu().detach().numpy() / count],
                                    "mean_token_accuracy": [mean_token_accuracy.cpu().detach().numpy() / count]})
        test_losses.to_csv(os.path.join(test_result_path, "test_losses.csv"), index=False)

        print(f"Testing Loss: {val_loss / count:.4f}")
        print(f"Average token reconstruction accuracy: {mean_token_accuracy / count:.4f}")

        # Plot predicted vs. true quantiles of 10 random functions
        number_of_samples = 10
        test_sample = test_pred_df.sample(number_of_samples, random_state=42)  # .transpose()
        melted_samples = []
        for col in test_sample.columns:
            # only melt quantile columns
            if test_sample[col].dtype != 'float32':
                continue
            col_sub = pd.DataFrame({"function": range(1, (number_of_samples + 1)),
                                    "quantile": float(col.split("Q")[1]),
                                    "parameter": test_sample[col]})
            # define split column
            if "pred" in col:
                col_sub["split"] = "prediction"
            else:
                col_sub["split"] = "true"
            melted_samples.append(col_sub)

        samples_plt_df = pd.concat(melted_samples).reset_index()
        palette = sns.color_palette("tab10", number_of_samples)
        sns.set_style("whitegrid")
        sns.set(rc={'figure.figsize': (14, 12)})
        sns.lineplot(data=samples_plt_df,
                     y="quantile",
                     x="parameter",
                     hue="function",
                     style="split",
                     style_order=["true", "prediction"],
                     palette=palette)
        plt.savefig(os.path.join(test_result_path, "quantile_sample_plot.png"))
        test_sample.to_csv(os.path.join(test_result_path, "quantile_sample_plot.csv"), index=False)

        # Plot in range 0-10
        sns.set_style("whitegrid")
        sns.set(rc={'figure.figsize': (14, 12)})
        sns.lineplot(data=samples_plt_df,
                     y="quantile",
                     x="parameter",
                     hue="function",
                     style="split",
                     style_order=["true", "prediction"],
                     palette=palette)
        plt.xlim(0, 10)
        plt.savefig(os.path.join(test_result_path, "quantile_sample_plot_range0-10.png"))
    ray.shutdown()