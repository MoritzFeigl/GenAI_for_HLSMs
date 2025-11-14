This repository provides the complete workflow accompanying the publication Distilling Hydrological and Land Surface Model Parameters from Physio-Geographical Properties Using Text-Generating AI. It includes data-preparation scripts, VAE and LSTM modeling pipelines, experiment orchestration tools, and analysis utilities. All scripts are designed for execution on a Linux HPC cluster with SLURM and may require adaptation (paths, resource settings, etc.) for other systems. A detailed list of required external raw datasets is provided below.

## Repository structure
* [code/](code)
  * [01_data_preparation/](code/01_data_preparation) scripts used for preparing all relevant data sets (spatial data and discharge time series) for mHM – requires raw data
  * [02_utility_functions/](code/02_utility_functions) GenAI_para utility functions used in the following steps
  * [03_experiment_setup/](code/03_experiment_setup) set up experiments including training/test split, run time estimation, and default parameter runs
  * [04_sensitivity_analysis/](code/04_sensitivity_analysis) sensitivity analysis of all optimizable mHM parameters (coefficients of mHM TFs) and consequent selection of relevant parameters for optimization
  * [05_VAE/](code/05_VAE)
    * [01_VAE_preprocessing/](code/05_VAE/01_VAE_preprocessing) scripts for generating and preparing data for VAE training
    * [02_VAE_training/](code/05_VAE/02_VAE_training) VAE scripts as the Python package 'GenAI_para' and scripts used for developing architecture, hyperparameter settings, and training of all VAEs (also containing a bash script for creating the correct conda environment)
  * [06_run_experiments/](code/06_run_experiments) setups and starting scripts for experiments
  * [07_analyze_results/](code/07_analyze_results) scripts for analyzing training/validation results and plotting parameter fields
  * [08_LSTM_modeling/](code/08_LSTM_modeling) LSTM data preparation, basin splitting, model training, and run configurations
  * [09_compare_data/](code/09_compare_data) scripts to prepare data for comparing parameter results (intercept and k_s)
* basins_results.csv – CSV file containing all modeling results from the publication
* study_basins.csv – CSV file containing all information on the used gauges and their data sources
## 1. System requirements

### Operating systems
- Linux HPC cluster with SLURM job scheduling. All submitted scripts assume a Linux filesystem layout and SLURM directives. This workflow has been tested on the Vienna Scientific Cluster. For local usage, use the provided list of libraries to re-create this setup.

### Software dependencies
- **R (\>= 4.0)** with the package set enumerated in the installation guide. These scripts are
  executed directly via `Rscript` or sourced within shell wrappers (`*.sh`) and represent the
  primary workflow engine for data preparation, sensitivity analysis, job submission, and
  post-processing.
- **mHM 5.11** (external Fortran code base) to run hydrological simulations; obtain from
  <https://git.ufz.de/mhm/mhm/-/releases/v5.11.0>. Data-preparation scripts expect mHM binaries and
  configuration templates to be available.
- **SLURM client tools** (`sbatch`, `squeue`, etc.) on the submit host.

The project also relies on two Conda environments tailored to specific stages of the workflow:

- **GenAI_VAE environment (Python 3.9)**【F:code/05_VAE/02_VAE_training/GenAI_para_VAEs/01_setup/00_conda_env_setup.sh†L1-L16】
  - PyTorch with CUDA 11.7 support (`pytorch`, `pytorch-cuda=11.7` via `pytorch`/`nvidia` channels)
  - Ray Tune/AIR/Data for hyperparameter search (`ray[tune]`, `ray[data]`, `ray[air]`)
  - Transformers, torchtext, hpbandster, ConfigSpace
  - Plotting and analysis utilities (`matplotlib`, `seaborn`, `plotly=5.11.0`)
  - Spreadsheet and text I/O (`openpyxl`, `chardet`, `charset-normalizer==2.1.0`)
  - General-purpose helpers (`packaging`)

- **LSTM environment (Python 3.10)**【F:code/08_LSTM_modeling/environment_cuda11_8.yml†L1-L28】
  - PyTorch with CUDA 11.8, pandas, xarray, scipy, numba, netCDF4, h5py
  - Developer tooling (pytest, pytest-cov, yapf, bokeh, jupyter)
  - Documentation extras (`sphinx`, `sphinx-rtd-theme`, `nbsphinx`, `nbsphinx-link`)
  - TensorBoard and scikit-learn (via pip)

### Tested configurations
- Python 3.9 + CUDA 11.7 for VAE experiments.
- Python 3.10 + CUDA 11.8 for LSTM baselines.
- SLURM-based HPC environment with NVIDIA V100/A100-class GPUs (minimum 16 GB GPU memory) for deep-learning workloads.
- Multi-core CPU nodes (\>= 16 vCPUs and 64 GB RAM) for data preparation and mHM optimization.
- 
### Required non-standard hardware
- NVIDIA GPU supporting CUDA 11.7/11.8 for VAE and LSTM training.
- High-throughput networked storage for large basin datasets and intermediate outputs.

### External data prerequisites
Most  scripts expect a shared `data/` directory (referenced as `data_dir` in R and Python code). Populate the following sub-directories with externally sourced datasets **before** running the pipeline. Paths that the repository itself generates (for example, VAE batches or trained checkpoints) are intentionally excluded from this list.

| Relative path under `data/` | Dataset description and acquisition guidance | Consumed by |
| --- | --- | --- |
| `required_data/discharge/` | Daily discharge observations per basin. Download from the Global Runoff Data Centre (GRDC, <https://www.bafg.de/GRDC>) or the relevant German state portal. Details on the gauges data sources are listed in [`study_basins.csv`](study_basins.csv). Some of these required direct requests to the operator and are not openly accessible. | mHM discharge preparation scripts in `code/01_data_preparation` (for example, `create_basin_discharge.py`). |
| `required_data/DWD_CDC_data/` | Monthly temperature and precipitation grids from the DWD Climate Data Center. Retrieve NetCDF tiles from <https://opendata.dwd.de/climate_environment/CDC/grids_germany/monthly/> and mirror the folder layout expected by `process_mat_map.py`. | Climate summary utilities in `code/01_data_preparation` (e.g., `process_mat_map.py`). |
| `required_data/forcings/` | Gridded forcing time series (tavg, PET, precipitation) for Germany supplied by the UFZ Department of Computational Hydrosystems (<https://www.ufz.de/index.php?en=41160>, ask for data up to 2019). | Forcing harmonisation and mHM input builders such as `create_basin_forcings.py`. |
| `required_data/forcings_UFZ/` | Template mHM forcing files from Feigl et al. (2022, <https://doi.org/10.1029/2022WR031966>) used to align naming conventions and metadata. | Forcing comparison workflows in `code/01_data_preparation` when mapping new forcings onto the UFZ templates. |
| `required_data/MOD15A3H_LAI_31468/` | Leaf Area Index (LAI) composites downloaded via Google Earth Engine (MOD15A3H collection). Export GeoTIFF tiles per basin and stage them according to the sub-folder structure referenced by `process_additional_modis_predictors.py`. | MODIS processing scripts in `code/01_data_preparation` (for example, `process_additional_modis_predictors.py`). |
| `required_data/static/` | Static spatial predictors per basin originating from Feigl et al. (2022, <https://doi.org/10.1029/2022WR031966>). | Spatial predictor assembly utilities such as `create_basin_mat.py` and `create_basin_lai.py`. |

If your deployment stores data outside the repository, adjust the `data_dir` substitutions in the provided scripts accordingly so that each script can locate its external inputs.
---

## 2. Installation guide

1. **Clone the repository**
   ```bash
   git clone https://github.com/<your-org>/GenAI_for_HLSMs.git
   cd GenAI_for_HLSMs
   ```

2. **Prepare data directories**
   - Mirror the directory layout referenced in the scripts (for example, `data/01_study_basins_info`, `data/03_discharge`, etc.).
   - Place raw shapefiles, forcing grids, and discharge time series into the matching folders on shared storage.

3. **Create the GenAI_VAE environment**
   ```bash
   conda create -n GenAI_VAE python=3.9 --yes
   conda activate GenAI_VAE
   conda install pytorch pytorch-cuda=11.7 -c pytorch -c nvidia --yes
   pip install "ray[tune]" "ray[data]" "ray[air]"
   conda install -c conda-forge packaging --yes
   pip install transformers hpbandster ConfigSpace chardet charset-normalizer==2.1.0
   conda install -c pytorch torchtext --yes
   conda install -c anaconda openpyxl seaborn --yes
   conda install -c conda-forge matplotlib --yes
   conda install -c plotly plotly=5.11.0 --yes
   ```

4. **Create the LSTM environment (optional)**
   ```bash
   conda env create -f code/08_LSTM_modeling/environment_cuda11_8.yml
   conda activate neuralhydrology
   ```

5. **Install R libraries**
   - Install the R packages referenced throughout the orchestration and analysis scripts:
     `RColorBrewer`, `RNetCDF`, `dplyr`, `ggplot2`, `gridExtra`, `magrittr`, `ncdf4`, `patchwork`,
     `raster`, `rasterVis`, `rnaturalearth`, `sensobol`, `sf`, `terra`, `viridis`, and `wesanderson`
     (plus optional helpers `ggthemes` and `stars` for specific plots). For example:
     ```r
     install.packages(c(
       "RColorBrewer", "RNetCDF", "dplyr", "ggplot2", "gridExtra", "magrittr",
       "ncdf4", "patchwork", "raster", "rasterVis", "rnaturalearth",
       "sensobol", "sf", "terra", "viridis", "wesanderson"
     ))
     ```
   - Configure SLURM submission commands (`sbatch`, `srun`) to be accessible from the login node.

A full setup on a typical research workstation (8-core CPU, 32 GB RAM, SSD storage, fast internet) takes approximately **30–45 minutes**, dominated by Conda environment creation and PyTorch downloads.

## 3. Demo

A lightweight sanity check is provided to validate your Python environment by summarizing the included `study_basins.csv` file.

1. Activate the GenAI_VAE (or any Python environment) and run:
   ```bash
   conda activate GenAI_VAE
   python - <<'PY'
   import csv
   import itertools
   from pathlib import Path
   path = Path('study_basins.csv')
   with path.open(encoding='latin-1') as f:
       reader = csv.DictReader(f)
       rows = list(itertools.islice(reader, 5))
   for row in rows:
       print({key: row[key] for key in ('Stat_ID','Station_Name','River','Start_Date','End_Date','split')})
   with path.open(encoding='latin-1') as f:
       count = sum(1 for _ in f) - 1
   print('\nNumber of basins:', count)
   PY
   ```

2. **Expected output**
   ```text
   {'Stat_ID': '6335115', 'Station_Name': 'GROLSHEIM', 'River': 'NAHE', 'Start_Date': '1972-11-01', 'End_Date': '2017-10-12', 'split': 'Validation'}
   {'Stat_ID': '6335125', 'Station_Name': 'SCHWAIBACH', 'River': 'KINZIG', 'Start_Date': '1914-01-01', 'End_Date': '2019-03-31', 'split': 'Training'}
   {'Stat_ID': '6335290', 'Station_Name': 'STEIN', 'River': 'KOCHER', 'Start_Date': '1911-06-01', 'End_Date': '2019-12-31', 'split': 'Validation'}
   {'Stat_ID': '6335304', 'Station_Name': 'FRANKFURT-AM-MAIN', 'River': 'MAIN', 'Start_Date': '1963-11-01', 'End_Date': '2019-12-31', 'split': 'Validation'}
   {'Stat_ID': '6335350', 'Station_Name': 'LEUN-(NEU)', 'River': 'LAHN', 'Start_Date': '1935-11-01', 'End_Date': '2019-12-31', 'split': 'Validation'}

   Number of basins: 162
   ```

3. **Run time**: < 10 seconds on a standard desktop (the script only parses a ~160-row CSV file).

---


## 4. Instructions for use

1. **Data preparation**
   - Configure paths in `code/01_data_preparation` scripts to point to your raw forcing, spatial, and discharge datasets.
   - Execute the R/py scripts to harmonize projections and generate mHM-ready inputs (e.g., run the utilities referenced in `code/02_utility_functions/04_first_run_preparation.sh`).

2. **VAE preprocessing and training**
   - Use `code/05_VAE/01_VAE_preprocessing` to create VAE training batches.
   - Launch dataset splitting with `01_setup/01_GenAI_vae_datasets.py`; adjust filesystem paths to your storage location before submitting through SLURM.
   - Submit training jobs defined in `code/05_VAE/02_VAE_training/GenAI_para_VAEs/05_ITCN_training/*.sh`, adapting SLURM resource requests (nodes, GPUs, walltime) to your cluster.

3. **Experiment orchestration**
   - Generate experiment directories and SLURM scripts with `code/06_run_experiments/01_start_experiments.R`. Update `main_path` to your workspace before sourcing the script in R; it will write run-specific submission scripts and call `sbatch` automatically.

4. **Hydrological model runs**
   - Ensure the mHM executable and configuration templates referenced in `code/02_utility_functions/05_run_mhm_preparation.sh` and related utilities point to your installation. These scripts synchronize transfer functions and trigger optimization runs through SLURM.

5. **Result analysis**
   - After experiments finish, use scripts in `code/07_analyze_results` to compute validation metrics and produce plots. Many scripts expect result directories under `results/` as created by the orchestration step.

6. **LSTM baseline (optional)**
   - Activate the `neuralhydrology` environment and follow the notebooks/scripts in `code/08_LSTM_modeling` to train comparison LSTM models. Adjust CUDA device IDs and data paths as necessary.

7. **Data comparison**
   - Run utilities in `code/09_compare_data` to prepare intercept and \(k_s\) parameter comparisons between GenAI-derived and reference datasets.

Adapt the provided SLURM scripts to match your scheduler’s partition names, account information, and filesystem layout. Because raw hydrological data are not bundled with the repository, you must supply all required shapefiles, forcing grids, and discharge records before executing the full workflow.

---

## Getting help
- Review inline comments within the SLURM submission scripts for resource recommendations.
- Open an issue in the repository with detailed logs if you encounter problems reproducing experiments.
