This repository contains the code of the publication *Distilling Hydrological and Land Surface Model Parameters from Physio-Geographical Properties Using Text-Generating AI*. A list of required raw data is given below. All scripts are written to work on a Linux cluster using SLURM and need to be adjusted (specifications, paths, etc.) for other clusters or local usage.

## Content of this repository
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

## Required data
* mHM model: the mesoscale Hydrological Model version 5.11 from <https://git.ufz.de/mhm/mhm/-/releases/v5.11.0>
* discharge: discharge data downloaded from GRDC or from different German state online services. Data for some gauges were only available after directly contacting the operator and are not publicly available. See [study_basins.csv](study_basins.csv) for details.
* DWD CDC data: monthly temperature and precipitation, downloadable from <https://opendata.dwd.de/climate_environment/CDC/grids_germany/monthly/>
* forcings: files with the forcings (tavg, pet, pre) for all of Germany provided by the UFZ Department of Computational Hydrosystems (<https://www.ufz.de/index.php?en=34211>)
* forcings_UFZ: mHM forcings from Feigl et al. (2022, <https://doi.org/10.1029/2022WR031966>) as templates
* MOD15A3H_LAI_31468: LAI data downloaded from Google Earth Engine
* static: original spatial predictors per basin from Feigl et al. (2022, <https://doi.org/10.1029/2022WR031966>)

