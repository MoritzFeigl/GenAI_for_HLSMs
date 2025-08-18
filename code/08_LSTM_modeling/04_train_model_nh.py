# Libararies
import os
import pickle
from pathlib import Path
import pandas as pd
from neuralhydrology.utils.config import Config
from neuralhydrology.nh_run import start_run, eval_run
from neuralhydrology.utils.nh_results_ensemble import create_results_ensemble
from neuralhydrology.evaluation.utils import metrics_to_dataframe
from neuralhydrology.evaluation.metrics import nse, kge

import numpy as np

main_path = Path("C:/Users/morit/Dropbox/Projekte/GenAI_para_for_HLSMs/code/08_LSTM_modeling")

# 1. initial run with all statics --------------------------------------------------------------------------------------
cfg = Config(Path("GenAI_para_ger_regional_training.yml"))
cfg = cfg.as_dict()
cfg["experiment_name"] = "all_statics"
main_run_dir = cfg["run_dir"]
cfg = Config(cfg)
cfg.dump_config(main_path, "all_statics.yml")
start_run(Path("all_statics.yml"))

# eval run
stored_runs = os.listdir(main_run_dir)
all_statics_run = [run for run in stored_runs if "all_statics" in run][0]
run_dir = Path(f"{main_run_dir}/{all_statics_run}")

# get best epoch
val_dirs = os.listdir(run_dir / "validation")
median_nse_vals = []
for dir in val_dirs:
    val_metrics = pd.read_csv(run_dir / "validation" / dir / "validation_metrics.csv")
    median_nse_vals.append(val_metrics["NSE"].median())
best_epoch = median_nse_vals.index(max(median_nse_vals)) + 1

# test evaluation
eval_run(run_dir=run_dir, period="test", epoch = best_epoch)



# 2. Ensemble runs with most important statics -------------------------------------------------------------------------
for i in range(8):
    print(f"run_{i+1}")
    if i == 0:
        continue
    cfg = Config(Path("study_regional_training.yml"))
    cfg = cfg.as_dict()
    cfg["static_attributes"] = ['map', 'aridity', 'area', 'dem', 'high_prec_dur', 'frac_snow']
    cfg["seed"] = cfg["seed"] + i + 1
    cfg["experiment_name"] = f"run_{i+1}"
    cfg = Config(cfg)

    cfg.dump_config(main_path, f"run_{i+1}.yml")
    start_run(Path(f"run_{i+1}.yml"))

    # eval run
    stored_runs = os.listdir(main_run_dir)
    current_run = [run for run in stored_runs if f"run_{i+1}" in run][0]
    run_dir = Path(f"{main_run_dir}/{current_run}")

    # get best epoch
    val_dirs = os.listdir(run_dir / "validation")
    median_nse_vals = []
    for dir in val_dirs:
        val_metrics = pd.read_csv(run_dir / "validation" / dir / "validation_metrics.csv")
        median_nse_vals.append(val_metrics["NSE"].median())
    best_epoch = median_nse_vals.index(max(median_nse_vals)) + 1

    # test evaluation
    eval_run(run_dir=run_dir, period="test", epoch=best_epoch)
    eval_run(run_dir=run_dir, period="train", epoch=best_epoch)
    eval_run(run_dir=run_dir, period="validation", epoch=best_epoch)


# create all test results
stored_runs = os.listdir(main_run_dir)
run_dirs = [Path(f"{main_run_dir}/" + run) for run in stored_runs if "run_" in run]


# Ensemble results
ensemble_results = create_results_ensemble(run_dirs, period='test')
ensemble_results_train = create_results_ensemble(run_dirs, period='train')
ensemble_results_val = create_results_ensemble(run_dirs, period='validation')

# helper function turn values to (m³/s)
def mm_to_m3_discharge(ser: pd.Series, area: float) -> pd.Series:
    """Helper function to reverse  discharge normalization by basin area"""
    return (ser / (1000 * 86400)) * (1e6 * area)

data_dir = Path("C:/Users/morit/Dropbox/Projekte/GenAI_para_for_HLSMs/data/09_lstm_data")
attributes = pd.read_csv(data_dir / "static" / "basin_attributes.csv", sep=";")
for basin in list(ensemble_results.keys()):
    basin_area = attributes["area"].loc[[str(basin_id) in basin for basin_id in attributes["basin"]]].item()
    ds = ensemble_results[basin]['1D']['xr']
    ds = ds.copy()
    for var in ['discharge_obs', 'discharge_sim']:
        ds[var] = mm_to_m3_discharge(ds[var], basin_area)
    ensemble_results[basin]['1D']['xr'] = ds

for basin in list(ensemble_results_train.keys()):
    basin_area = attributes["area"].loc[[str(basin_id) in basin for basin_id in attributes["basin"]]].item()
    ds = ensemble_results_train[basin]['1D']['xr']
    ds = ds.copy()
    for var in ['discharge_obs', 'discharge_sim']:
        ds[var] = mm_to_m3_discharge(ds[var], basin_area)
    ensemble_results_train[basin]['1D']['xr'] = ds

for basin in list(ensemble_results_val.keys()):
    basin_area = attributes["area"].loc[[str(basin_id) in basin for basin_id in attributes["basin"]]].item()
    ds = ensemble_results_val[basin]['1D']['xr']
    ds = ds.copy()
    for var in ['discharge_obs', 'discharge_sim']:
        ds[var] = mm_to_m3_discharge(ds[var], basin_area)
    ensemble_results_val[basin]['1D']['xr'] = ds

def update_nse_kge_for_basin(results: dict, basin: str, model: str = '1D'):
    """
    Recompute NSE and KGE for one basin/model,
    based on whatever discharge_obs/discharge_sim are currently in the xr‑Dataset.
    """
    # grab the little dict and its Dataset
    sub = results[basin][model]
    ds  = sub['xr']

    # extract the two series (as pandas Series, since your nse/kge expect pd.Series)
    obs = ds['discharge_obs'].to_series()
    sim = ds['discharge_sim'].to_series()

    # call your metrics
    new_nse = nse(obs, sim)
    new_kge = kge(obs, sim)

    # overwrite
    sub['NSE'] = new_nse
    sub['KGE'] = new_kge

for b in ensemble_results:
    update_nse_kge_for_basin(ensemble_results, b, model='1D')

for b in ensemble_results_train:
    update_nse_kge_for_basin(ensemble_results_train, b, model='1D')
for b in ensemble_results_val:
    update_nse_kge_for_basin(ensemble_results_val, b, model='1D')

# store results
output_dir = main_run_dir / "ensemble_run"
os.makedirs(output_dir, exist_ok=True)
file_name = output_dir / f"test_ensemble_results.p"
pickle.dump(ensemble_results, open(file_name, 'wb'))


# Evaluate ensemble results: NSE and KGE
cfg = Config(run_dirs[0] / 'config.yml')
df_nse_kge = metrics_to_dataframe(ensemble_results, ["NSE", "KGE"], cfg.target_variables)
df_nse_kge_train = metrics_to_dataframe(ensemble_results_train, ["NSE", "KGE"], cfg.target_variables)
df_nse_kge_val = metrics_to_dataframe(ensemble_results_val, ["NSE", "KGE"], cfg.target_variables)

full_df = pd.concat([df_nse_kge, df_nse_kge_train, df_nse_kge_val], axis=0)
full_df.to_csv(output_dir / f"full_ensemble_metrics.csv")
# Evaluate ensemble results: log-NSE

# compute log discharge
for basin in list(ensemble_results.keys()):
    ds = ensemble_results[basin]['1D']['xr']
    ds = ds.copy()
    for var in ['discharge_obs', 'discharge_sim']:
        ds[var] = np.log(ds[var])
    ensemble_results[basin]['1D']['xr'] = ds

for b in ensemble_results:
    update_nse_kge_for_basin(ensemble_results, b, model='1D')


# Evaluate log NSE and concat metrics
df_lnse = metrics_to_dataframe(ensemble_results, ["NSE"], cfg.target_variables)
df_lnse.columns = ["lNSE"]
df = pd.concat([df_nse_kge, df_lnse], axis=1)

# store results
file_name = output_dir / f"test_ensemble_metrics.csv"
df.to_csv(file_name)
