# mHM ET sensitivity analysis
# GenAI_para-ET-mHM project

# define run name to create environment automatically
run_name <- "04_SA"

# lib path and libraries
library(ggplot2)
library(sensobol)

# create run and result folders
run_dir <- paste0("/gpfs/data/fs71468/GenAI_para_runs/runs/", run_name)
result_dir <- paste0("/gpfs/data/fs71468/GenAI_para_runs/results/", run_name)
code_dir <- "/gpfs/data/fs71468/GenAI_para_runs/code/"
data_dir <- "/gpfs/data/fs71468/GenAI_para_runs/data/"
dir.create("/gpfs/data/fs71468/GenAI_para_runs/runs", showWarnings = FALSE)
dir.create(run_dir, showWarnings = FALSE)
dir.create("/gpfs/data/fs71468/GenAI_para_runs/results", showWarnings = FALSE)
dir.create(result_dir, showWarnings = FALSE)
setwd(run_dir)
source(paste0(code_dir, "02_utility_functions/01_GenAI_para_utility_functions.R"))


# 2. Sobol samples for 59 Parameters --------------------------------------------------------
# get parameter infos
parameters <- read_nml(paste0(data_dir, "/05_mhm/master_mhm_config/mhm_parameter.nml"))
parameter_df <- do.call(rbind, parameters$mhm_parameters)
ind_of_rel_paras <- c(1:58, 68:90)
num_para_bounds <- parameter_df[ind_of_rel_paras, c(1, 2)]
standard_parameters <- parameter_df[ind_of_rel_paras, 3]

# Create sample matrix to compute first and total-order indices:
if(!file.exists(paste0(result_dir, "/01_initial_sobol_samples.rds"))){
  params <- names(standard_parameters)
  N <- 40
  mat <- sobol_matrices(N = N, params = params, type = "QRN")
  cat("estimated time: ", nrow(mat)*0.057/60/24, " days or ", nrow(mat)*0.057/60, " hours")
  for(ipara in 1:ncol(mat)){
    mat[, ipara] <- qunif(mat[, ipara], num_para_bounds[ipara, 1], num_para_bounds[ipara, 2])
  }
  saveRDS(mat, paste0(result_dir, "/01_initial_sobol_samples.rds"))
}

# 2. Inital SA with 59 Parameters --------------------------------------------------------

# create SA run environment
if(!dir.exists(paste0(run_dir, "/mhm"))){
  
  # chooose basins for sensitivity analysis
  # one of each N/S/E/W and central Germany
  # smaller basins are preferred to reduce run time
  
  basins <- read.csv(paste0(data_dir, "06_training_validation_data/basin_split.csv"))
  basins <- basins[basins$split == "Training",]
  extent <- data.frame(Y = c(min(basins$Norting_Y), max(basins$Norting_Y)),
                       X = c(min(basins$Easting_X), max(basins$Easting_X)))
  # split into 3x3 grid and choose the smallest basin in relevant cells
  y_step <- (extent$Y[2] - extent$Y[1])/3
  x_step <- (extent$X[2] - extent$X[1])/3
  
  # north west
  N_basins <- basins[basins$Easting_X > extent$X[1] + x_step  & 
                       basins$Easting_X < extent$X[2] - x_step &
                       basins$Norting_Y > extent$Y[1] + 2*y_step & 
                       basins$Norting_Y < extent$Y[2], ]
  N_basin <- N_basins[N_basins$area == min(N_basins$area), ]
  # South
  S_basins <- basins[basins$Easting_X > extent$X[1] + x_step & 
                       basins$Easting_X < extent$X[2] - x_step &
                       basins$Norting_Y > extent$Y[1] & 
                       basins$Norting_Y < extent$Y[2] - 2*y_step, ]
  S_basin <- S_basins[S_basins$area == min(S_basins$area), ]
  # East
  E_basins <- basins[basins$Easting_X > extent$X[1] + 2*x_step & 
                       basins$Easting_X < extent$X[2] &
                       basins$Norting_Y > extent$Y[1] + y_step & 
                       basins$Norting_Y < extent$Y[2] - y_step, ]
  E_basin <- E_basins[E_basins$area == min(E_basins$area), ]
  # West
  W_basins <- basins[basins$Easting_X > extent$X[1] & 
                       basins$Easting_X < extent$X[2] - 2*x_step &
                       basins$Norting_Y > extent$Y[1] + y_step & 
                       basins$Norting_Y < extent$Y[2] - y_step, ]
  W_basin <- W_basins[W_basins$area == min(W_basins$area), ]
  # Central
  C_basins <- basins[basins$Easting_X > extent$X[1] + x_step & 
                       basins$Easting_X < extent$X[2] - x_step &
                       basins$Norting_Y > extent$Y[1] + y_step & 
                       basins$Norting_Y < extent$Y[2] - y_step, ]
  C_basin <- C_basins[C_basins$area == min(C_basins$area), ]
  selected_basins <- cbind(rbind(N_basin, S_basin, E_basin, W_basin, C_basin),
                           selected = TRUE)
  # Plot selected basins
  ggplot(rbind(cbind(basins, selected = FALSE), selected_basins),
         aes(Easting_X, Norting_Y, color = selected)) + geom_point() +
    geom_hline(yintercept = extent$Y[1] + y_step) + 
    geom_hline(yintercept = extent$Y[2] - y_step)  +
    geom_vline(xintercept = extent$X[1] + x_step) + 
    geom_vline(xintercept = extent$X[2] - x_step) 
  ggsave(paste0(result_dir, "/02_SA_basins_position.png"))
  # Ems
  # FICHTENBERGER-ROT -> Kocher
  # MUEGLITZ -> Elbe
  # SCHWARZBACH -> Blies
  # Gersprenz -> Main
  write.csv(selected_basins, paste0(result_dir, "/03_selected_basins.csv"), 
            row.names = FALSE)
  
  # mhm
  file.copy(paste0(data_dir, "/05_mhm/mhm"), run_dir, 
            recursive=TRUE, overwrite = TRUE)
  # remove build folder
  try(unlink(paste0(run_dir, "/mhm/build"), recursive = TRUE))
  
  # copy and adapt first_run_preparation script & run_mhm_preparation script
  # first run compilation script
  run_preparation <- readLines(paste0(code_dir, "02_utility_functions/04_first_run_preparation.sh"))
  run_preparation <- gsub("PROJECTPATH=/home/fs71468/mfeigl/GenAI_para_for_HLSMs/runs/",
                          paste0("PROJECTPATH=", run_dir), run_preparation, fixed = TRUE)
  run_preparation <- gsub("6335125", selected_basins$Stat_ID[1], run_preparation)
  run_preparation <- gsub("<DATA_DIR>", data_dir, run_preparation, fixed = TRUE)
  writeLines(run_preparation, paste0(run_dir, "/04_first_run_preparation.sh"))
  
  # run preparation (mhm compilation) script
  run_preparation <- readLines(paste0(code_dir, "02_utility_functions/05_run_mhm_preparation.sh"))
  run_preparation <- gsub("PROJECTPATH=/home/fs71468/mfeigl/GenAI_para_for_HLSMs/runs/",
                          paste0("PROJECTPATH=", run_dir), run_preparation, fixed = TRUE)
  run_preparation <- gsub("6335125", selected_basins$Stat_ID[1], run_preparation, fixed = TRUE)
  run_preparation <- gsub("<DATA_DIR>", data_dir, run_preparation, fixed = TRUE)
  writeLines(run_preparation, paste0(run_dir, "/05_run_mhm_preparation.sh"))
  
  # create output directories
  dir.create(paste0(run_dir, "/Training"), showWarnings = FALSE)
  
  # adapt paths in mhm.nml for each basin and change time period
  for(basin in selected_basins$Stat_ID){
    # create basin specific mhm config and output folder
    dir.create(paste0(run_dir, "/Training/", basin), showWarnings = FALSE)
    basin_dir <- paste0(run_dir, "/Training/", basin)
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
    mpr_nml <- sapply(mpr_nml, function(x) gsub("<DATA_DIR>", data_dir, x))
    mpr_nml <- sapply(mpr_nml, function(x) gsub("<SPLIT>", "Training", x))
    mpr_nml <- sapply(mpr_nml, function(x) gsub("<BASIN>", basin, x))
    writeLines(mpr_nml, paste0(basin_dir,"/config/mpr.nml"))
    
    # create basin output directory
    file.copy(paste0(data_dir, "06_training_validation_data/Training/", basin, "/output"),
              basin_dir, recursive=TRUE, overwrite = TRUE)
    
  }
  
  # Change output to only aET monthly
  state_and_fluxes[["aET"]] <- TRUE
  change_mhm_outputs(basins_dir="Training", state_and_fluxes, time_step = "monthly")
  
  # Create run scrpts
  dir.create("training_run_scripts", showWarnings = FALSE)
  for(basin in selected_basins$Stat_ID){
    run_script_name <- paste0(run_dir, "/training_run_scripts/basin_" ,basin, ".sh")
    fileConn <- file(run_script_name)
    writeLines(
      c("#!/bin/bash",
        paste0("cd ", run_dir, "/Training/", basin, "/config"),
        "./mhm > output.txt  2>&1"),
      fileConn)
    close(fileConn)
  }
  
  # Set standard mHM parameter
  KSat <- "PTF_Ks_curveSlope * exp((PTF_Ks_constant + PTF_Ks_sand * sand - PTF_Ks_clay * clay) * log(Ks_c_base))"
  FieldCap <- "ThetaS * exp(FieldCap_c1 * (FieldCap_c2 + log10(KSat)) * log(vGenu_n))"
  froots_1 <- "(1.0 - rootFractionCoefficient_forest ** (z_upper_bound * 100.0)) - (1.0 - rootFractionCoefficient_forest ** (z_lower_bound * 100.0))"
  froots_2 <- "(1.0 - (rootFractionCoefficient_forest - rootFractionCoefficient_pervious)) ** (z_upper_bound * 100.0) - (1.0 - (rootFractionCoefficient_forest - rootFractionCoefficient_pervious) ** (z_lower_bound * 100.0))"
  ThetaS_1 <- "PTF_lower66_5_constant + PTF_lower66_5_clay * clay + PTF_lower66_5_Db * bd"
  ThetaS_2 <- "PTF_higher66_5_constant + PTF_higher66_5_clay * clay + PTF_higher66_5_Db * bd"
  Canopy_Intercept <- "canopyInterceptionFactor * lai_class"
  
  parameter_list <- list("Ksat" = KSat,
                         "FieldCap" = FieldCap,
                         "fRoots_1" = froots_1, "fRoots_2" = froots_2,
                         "ThetaS_1" = ThetaS_1, "ThetaS_2" = ThetaS_2,
                         "L1_Max_Canopy_Intercept" = Canopy_Intercept)
  
  # update TFs and numeric TF coefficients
  update_tfs_in_mpr_nml(run_dir, setup=NULL, parameter_list=parameter_list, scaling = FALSE)
  
  # # initial compilation
  system("bash 04_first_run_preparation.sh > 04_first_run_preparation.out & wait")
  
}


# which numeric parameter to change in SA
parameter_nml <- read_nml(paste0(data_dir, "/05_mhm/master_mhm_config/mhm_parameter.nml"))
parameter_df <- do.call(rbind, parameter_nml$mhm_parameters)
ind_of_rel_paras <- c(1:58, 68:90)
num_para_bounds <- parameter_df[ind_of_rel_paras, c(1, 2)]
standard_parameters <- parameter_df[ind_of_rel_paras, 3]

# define objective function for SA
objective_function <- function(numeric_parameters){
  
  # write new mpr.nml and mhm_parameter.nml for all basins
  update_numeric_parameters(numeric_parameters=numeric_parameters, setup = NULL,
                            run_dir=run_dir,
                            parameter_nml=parameter_nml, 
                            ind_of_num_paras=ind_of_rel_paras)
  
  # run mhm with numeric parameters
  start_eval <- Sys.time()
  system(paste0("export LD_LIBRARY_PATH=$LIBRARY_PATH
  cd training_run_scripts", "
                  for j in *
                  do
                  bash $j &
                  done
                  wait"))
  end_eval <- Sys.time()
  cat("Evaluation took:\n")
  print(end_eval - start_eval)
  
  # get ET and Q results and compute quality criteria
  basin_qc <- basins_quality_criteria(run_dir, setup = NULL, data_dir)
  loss_list <- compute_loss(basin_qc, mode="QET")
  loss <- loss_list[[1]]
  domain_wloss <- loss_list[[2]]
  
  # write mean criteria values to csv
  results <- apply(basin_qc[, -1], 2, mean)
  results <- as.data.frame(matrix(c(results, loss), ncol = 5))
  results <- cbind(results, "time" = Sys.time())
  names(results) <- c(colnames(basin_qc)[-1], "loss", "time")
  results <- cbind(results, matrix(numeric_parameters, nrow = 1))
  names(results)[7:ncol(results)] <- names(numeric_parameters)
  
  if("04_initial_SA_results.csv" %in% list.files(result_dir)){
    all_results <- read.csv(paste0(result_dir, "/04_initial_SA_results.csv"))
    names(results) <- names(all_results)
    all_results <- rbind(all_results, results)
  } else {
    all_results <- results
  }
  write.csv(all_results, paste0(result_dir, "/04_initial_SA_results.csv"), row.names = FALSE)
  return(loss)
}

# define eval function
run_objfun <- function(dt){
  return(apply(dt, 1, objective_function))
}

# load precomputed sobol samples
mat <- readRDS(paste0(result_dir, "/01_initial_sobol_samples.rds"))

# remove any previous results
try({file.remove(paste0(result_dir, "/04_initial_SA_results.csv"))})

# compute sobol samples results
start_run <- Sys.time()
Y <- run_objfun(mat)
end_run <- Sys.time()
cat("Full Run took:\n")
print(end_run - start_run)

# store results
saveRDS(Y, paste0(result_dir, "/05_initial_sobol_output.rds"))

# start follow up sensitivity analysis with reduced paramerter set
setwd(paste0(code_dir, "04_sensitivity_analysis"))
system("sbatch 02_SA.sh")
