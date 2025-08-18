
library(magrittr)
library(ggplot2)

# Create Training environment ----------------------------------------------------------
run_name <- "Training"
main_path <- "/gpfs/data/fs71468/GenAI_para_runs/"
run_dir <- paste0(main_path, "runs/", run_name)
result_dir <- paste0(main_path, "results/", run_name)
code_dir <- paste0(main_path, "code/")
data_dir <- paste0(main_path, "data/")
setup <- n_parallel_setups <- 1
dir.create(run_dir, showWarnings = FALSE)
dir.create(result_dir, showWarnings = FALSE)
setwd(run_dir)
source(paste0(code_dir, "02_utility_functions/01_GenAI_para_utility_functions.R"))
# basin list
basins <- list.files(paste0(data_dir, "06_training_validation_data/Training"))
# set up mHM environment
if(!dir.exists(paste0(run_dir, "/mhm_", setup))){
  # mhm
  try(unlink(paste0(run_dir, "/mhm_", setup), recursive = TRUE))
  file.copy(paste0(data_dir, "/05_mhm/mhm"), run_dir, 
            recursive = TRUE, overwrite = TRUE)
  file.rename(paste0(run_dir, "/mhm"), paste0(run_dir, "/mhm_", setup))
  try(unlink(paste0(run_dir, "/mhm_", setup, "/build"), recursive = TRUE))
  
  # first run compilation script
  run_preparation <- readLines(paste0(code_dir, "02_utility_functions/04_first_run_preparation.sh"))
  run_preparation <- gsub("PROJECTPATH=/home/fs71468/mfeigl/GenAI_para_for_HLSMs/runs/",
                          paste0("PROJECTPATH=", run_dir), run_preparation, fixed = TRUE)
  run_preparation <- gsub("/mhm/", paste0("/mhm_", setup, "/"), run_preparation, fixed = TRUE)
  run_preparation <- gsub("ROJECTPATH/Training", paste0("ROJECTPATH/Training_", setup), run_preparation, fixed = TRUE)
  run_preparation <- gsub("6335125", basins[1], run_preparation)
  run_preparation <- gsub("<DATA_DIR>", data_dir, run_preparation, fixed = TRUE)
  writeLines(run_preparation, paste0(run_dir, "/04_first_run_preparation_", setup, ".sh"))
  
  # run preparation (mhm compilation) script
  run_preparation <- readLines(paste0(code_dir, "02_utility_functions/05_run_mhm_preparation.sh"))
  run_preparation <- gsub("PROJECTPATH=/home/fs71468/mfeigl/GenAI_para_for_HLSMs/runs/",
                          paste0("PROJECTPATH=", run_dir), run_preparation, fixed = TRUE)
  run_preparation <- gsub("/mhm/", paste0("/mhm_", setup, "/"), run_preparation, fixed = TRUE)
  run_preparation <- gsub("ROJECTPATH/Training", paste0("ROJECTPATH/Training_", setup), run_preparation, fixed = TRUE)
  run_preparation <- gsub("6335125", basins[1], run_preparation, fixed = TRUE)
  run_preparation <- gsub("<DATA_DIR>", data_dir, run_preparation, fixed = TRUE)
  writeLines(run_preparation, paste0(run_dir, "/05_run_mhm_preparation_", setup, ".sh"))
  
  # create basin mhm dirs
  dir.create(paste0(run_dir, "/Training_", setup), showWarnings = FALSE)
  
  # adapt paths in mhm.nml for each basin and change time period
  for(basin in basins){
    # create basin specific mhm config and output folder
    basin_dir <- paste0(run_dir, "/Training_", setup, "/", basin)
    dir.create(basin_dir, showWarnings = FALSE)
    file.copy(paste0(data_dir, "06_training_validation_data/Training/", basin, "/config/"),
              basin_dir, recursive=TRUE, overwrite = TRUE)
    
    # mhm.nml
    mhm_nml <- readLines(paste0(basin_dir,"/config/mhm.nml"))
    mhm_nml <- sapply(mhm_nml, function(x) gsub("<DATA_DIR>", data_dir, x))
    mhm_nml <- sapply(mhm_nml, function(x) gsub("<SPLIT>", "Training", x))
    mhm_nml <- sapply(mhm_nml, function(x) gsub("<BASIN>", basin, x))
    writeLines(mhm_nml, paste0(basin_dir,"/config/mhm.nml"))
    
    # mpr.nml
    file.copy(paste0(data_dir, "05_mhm/master_mhm_config/", "mpr.nml"),
              paste0(basin_dir,"/config/mpr.nml"), overwrite = TRUE)
    mpr_nml <- readLines(paste0(basin_dir,"/config/mpr.nml"))
    mpr_nml <- sapply(mpr_nml, function(x) gsub("<DATA_DIR>", data_dir, x,))
    mpr_nml <- sapply(mpr_nml, function(x) gsub("<SPLIT>", "Training", x))
    mpr_nml <- sapply(mpr_nml, function(x) gsub("<BASIN>", basin, x))
    writeLines(mpr_nml, paste0(basin_dir,"/config/mpr.nml"))
    
    # create basin output directory
    file.copy(paste0(data_dir, "06_training_validation_data/Training/", basin, "/output"),
              basin_dir, recursive=TRUE, overwrite = TRUE)
    
  }
  
  # Create run scripts
  dir.create(paste0("training_run_scripts_", setup), showWarnings = FALSE)
  for(basin in list.files(paste0(run_dir, "/Training_", setup))){
    run_script_name <- paste0(run_dir, "/training_run_scripts_", setup, 
                              "/basin_", basin, ".sh")
    fileConn <- file(run_script_name)
    writeLines(
      c("#!/bin/bash",
        paste0("cd ", run_dir, "/Training_", setup, "/", basin, "/config"),
        "./mhm > output.txt  2>&1"),
      fileConn)
    close(fileConn)
  }
  
  # # initial compilation
  comp_cmd <- paste0(paste0(
    "bash 04_first_run_preparation_",
    1:n_parallel_setups, ".sh & ", collapse = ""), "wait")
  
  system(comp_cmd)
}

# rerun all best Training results --------------------------------------------------------
for(experiment in c(0, 1, 2)){
  training_results <- read.csv("/gpfs/data/fs71468/GenAI_para_runs/results/training_results_2024-01-16.csv")
  exp_results <- unique(training_results[training_results$experiment == experiment, c("run", "mean_loss", "time")])
  run <- exp_results[exp_results$mean_loss == max(exp_results$mean_loss), "run"]
  
  
  # define experiment configuration
  source(paste0(code_dir, "06_run_experiments/experiment_configs/experiment_", experiment, "_setup.R"))
  
  # define mhm output --> always aET for result analysis
  state_and_fluxes[["aET"]] <- TRUE
  for(setup in 1:n_parallel_setups){
    change_mhm_outputs(basins_dir = paste0(run_dir, "/Training_", setup), 
                       state_and_fluxes, 
                       time_step = "monthly")
  }
  
  # get best parameter set
  path <- paste0(result_dir, "/../exp", experiment, "-run", run, "/")
  tf1 <- read.csv(paste0(path, "best_TFs_tracker_1.csv"))
  tf2 <- read.csv(paste0(path, "best_TFs_tracker_2.csv"))
  last_tfs <- rbind(tf1[nrow(tf1), ], tf2[nrow(tf2), ])
  best_parameter <- last_tfs[which.max(c(tail(tf1$loss, 1), tail(tf2$loss, 1))), ]
  
  
  if(experiment == 0){
    # Standard parameter
    KSat <- "PTF_Ks_curveSlope * exp((PTF_Ks_constant + PTF_Ks_sand * sand - PTF_Ks_clay * clay) * log(Ks_c_base))"
    FieldCap <- "ThetaS * exp(FieldCap_c1 * (FieldCap_c2 + log10(KSat)) * log(vGenu_n))"
    fRoots_1 <- "(1.0 - rootFractionCoefficient_forest ** (z_upper_bound * 100.0)) - (1.0 - rootFractionCoefficient_forest ** (z_lower_bound * 100.0))"
    fRoots_2 <- "(1.0 - (rootFractionCoefficient_forest - rootFractionCoefficient_pervious)) ** (z_upper_bound * 100.0) - (1.0 - (rootFractionCoefficient_forest - rootFractionCoefficient_pervious) ** (z_lower_bound * 100.0))"
    ThetaS_1 <- "PTF_lower66_5_constant + PTF_lower66_5_clay * clay + PTF_lower66_5_Db * bd"
    ThetaS_2 <- "PTF_higher66_5_constant + PTF_higher66_5_clay * clay + PTF_higher66_5_Db * bd"
    Canopy_Intercept <- "canopyInterceptionFactor * lai_class"
  } else {
    KSat <- best_parameter$Ksat
    FieldCap <- best_parameter$FieldCap
    fRoots_1 <- best_parameter$fRoots_1
    fRoots_2 <- best_parameter$fRoots_2
    ThetaS_1 <- best_parameter$ThetaS_1
    ThetaS_2 <- best_parameter$ThetaS_2
    Canopy_Intercept <- best_parameter$L1_Max_Canopy_Intercept
  }
  parameter_list <- list("Ksat" = KSat,
                         "FieldCap" = FieldCap,
                         "fRoots_1" = fRoots_1, "fRoots_2" = fRoots_2,
                         "ThetaS_1" = ThetaS_1, "ThetaS_2" = ThetaS_2,
                         "L1_Max_Canopy_Intercept" = Canopy_Intercept)
  
  
  
  # if all TFs are valid run mHM and compute losses
  trying <- try({update_tfs_in_mpr_nml(run_dir, setup, parameter_list, scaling = FALSE, 
                                       save_para_as_ncdf = TRUE, Training=TRUE)})
  if(class(trying) == "try-error"){
    for(basin in basins){
      basin_dir <- paste0(run_dir, "/Training_", setup, "/", basin)
      file.copy(paste0(main_path, "runs/exp", experiment, "-run", run, "/Training_2/", basin, "/config/mpr.nml"),
                paste0(basin_dir,"/config/mpr.nml"), overwrite = TRUE)
    }
    update_tfs_in_mpr_nml(run_dir, setup, parameter_list, scaling = FALSE, 
                          save_para_as_ncdf = TRUE, Training=TRUE)
  }
  
  # select numeric mhm para based on best results
  name_fix_df <- data.frame(matrix(NA, nrow=1, ncol=length(names(standard_parameters))))
  names(name_fix_df) <- names(standard_parameters)
  relevant_parameters <- names(data.frame(name_fix_df))
  best_num_para <- best_parameter[, relevant_parameters]
  
  # numeric para update
  update_numeric_parameters(best_num_para, run_dir, setup, 
                            parameter_nml, ind_of_num_paras, Training=TRUE)
  # Run mHM
  run_mhm(run_dir, setup, mhm_preparation = TRUE, Training=TRUE)
  
  # compute quality criteria and loss
  basin_qc <- basins_quality_criteria(run_dir, setup, data_dir, Training=TRUE, ET=TRUE)
  
  # save results, parameter fields and fluxes
  exp_result_dir <- paste0(main_path, "results/", run_name, "/exp", experiment, "-run", run, "/")
  dir.create(exp_result_dir)
  
  # 1 save basin_qc and loss
  write.csv(basin_qc, paste0(exp_result_dir, "basins_quality_criteria_training.csv"))
  # 2. save parameter files
  for(basin in basins){
    system(paste0("cp ", 
                  paste0(run_dir, "/Training_1/", basin, "/output/mHM_parameters.nc "),
                  paste0(exp_result_dir, basin, "_parameters.nc")))
    
    file.copy(paste0(run_dir, "/Training_1/", basin, "/output/mHM_parameters.nc"),
              paste0(exp_result_dir, basin, "_parameters.nc"), overwrite = TRUE)
  }
  # 3. save fluxes
  for(basin in basins){
    file.copy(paste0(run_dir, "/Training_1/", basin, "/output/mHM_Fluxes_States.nc"),
              paste0(exp_result_dir, basin, "_Fluxes_States.nc"), overwrite = TRUE)
  }
}



















