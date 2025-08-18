# Sensitivity analysis evaluation
# GenAI_para-ET-mHM project

# define run name to create environment automatically
run_name <- "04_SA"


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
library(sensobol)

# 1. Prepare sobol sample with reduced parameter set -------------------------------------
# get results from initial SA
mat <- readRDS(paste0(result_dir, "/01_initial_sobol_samples.rds"))
y <- readRDS(paste0(result_dir, "/05_initial_sobol_output.rds"))

# get parameter infos
parameters <- read_nml(paste0(data_dir, "/05_mhm/master_mhm_config/mhm_parameter.nml"))
parameter_df <- do.call(rbind, parameters$mhm_parameters)
ind_of_rel_paras <- c(1:58, 68:90)
num_para_bounds <- parameter_df[ind_of_rel_paras, c(1, 2)]
standard_parameters <- parameter_df[ind_of_rel_paras, 3]


# Compute and bootstrap the Sobol' indices:
R <- 100
params <- names(standard_parameters)
N <- 40
ind <- sobol_indices(Y = y, N = N, params = params,
                     first = "saltelli", total = "jansen", 
                     boot = TRUE, R = R, 
                     conf = 0.95, type = "norm"
)
# save computed indices
write.csv(ind$results, paste0(result_dir, "/06_initial_sobol_indices.csv"), row.names = FALSE)

# 2. SA with only sensitive Parameters ------------------------------------------------

# get initial SA results
results <- read.csv(paste0(result_dir, "/06_initial_sobol_indices.csv"))

# reduce parameters using the total sobol index
ti <- results[results$sensitivity == "Ti", c("original", "std.error", "parameters")]
# remove parameters with sensitivity < 10e-5
ti_reduced <- ti[ti$original > 10e-5, ]

# reduced parameter and parameter bounds
red_standard_parameters <- standard_parameters[names(standard_parameters) %in% ti_reduced$parameters]
red_num_para_bounds <- num_para_bounds[rownames(num_para_bounds) %in% ti_reduced$parameters, ]
saveRDS(red_standard_parameters, paste0(result_dir, "/07_red_standard_parameters.rds"))
saveRDS(red_num_para_bounds, paste0(result_dir, "/07_red_num_para_bounds.rds"))

# Create sample matrix to compute first and total-order indices:
red_params <- names(red_standard_parameters)
N <-  100
red_mat <- sobol_matrices(N = N, params = red_params, type = "QRN")
cat("estimated time: ", nrow(red_mat)*0.38/60/24, " days or ", nrow(red_mat)*0.38/60, " hours")
for(ipara in 1:ncol(red_mat)){
  red_mat[, ipara] <- qunif(red_mat[, ipara], 
                            red_num_para_bounds[ipara, 1], 
                            red_num_para_bounds[ipara, 2])
}
saveRDS(red_mat, paste0(result_dir, "/08_sobol_samples.rds"))

# 2. run SA ------------------------------------------------------------------------------

# which numeric parameter to change in SA
parameter_nml <- read_nml(paste0(data_dir, "/05_mhm/master_mhm_config/mhm_parameter.nml"))
parameter_df <- do.call(rbind, parameter_nml$mhm_parameters)

# run with reduced parameter set
standard_parameters <- readRDS(paste0(result_dir, "/07_red_standard_parameters.rds"))
num_para_bounds <- readRDS(paste0(result_dir, "/07_red_num_para_bounds.rds"))
ind_of_rel_paras <- which(rownames(parameter_df) %in% names(standard_parameters))

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
  if("09_SA_results.csv" %in% list.files(result_dir)){
    all_results <- read.csv(paste0(result_dir, "/09_SA_results.csv"))
    names(results) <- names(all_results)
    all_results <- rbind(all_results, results)
  } else {
    all_results <- results
  }
  write.csv(all_results, paste0(result_dir, "/09_SA_results.csv"), row.names = FALSE)
  return(loss)
}


# define eval function
run_objfun <- function(dt){
  return(apply(dt, 1, objective_function))
}

# load predcomputed sobol samples
mat <- readRDS(paste0(result_dir, "/08_sobol_samples.rds"))

# remove any previous results
try({file.remove(paste0(result_dir, "/09_SA_results.csv"))})

# compute sobol samples results
start_run <- Sys.time()
Y <- run_objfun(mat)
end_run <- Sys.time()
cat("Full Run took:\n")
print(end_run - start_run)

# store results
saveRDS(Y, paste0(result_dir, "/10_sobol_output.rds"))

# Compute sobol indices ------------------------------------------------------------------
red_mat <- readRDS(paste0(result_dir, "/08_sobol_samples.rds"))
y <- readRDS(paste0(result_dir, "/10_sobol_output.rds"))
red_standard_parameters <- readRDS(paste0(result_dir, "/07_red_standard_parameters.rds"))
red_params <- names(red_standard_parameters)


# Compute and bootstrap the Sobol' indices:
N <-  100
R <- 100
ind <- sobol_indices(Y = y, N = N, params = red_params,
                     first = "saltelli", total = "jansen", 
                     boot = TRUE, R = R, 
                     conf = 0.95, type = "norm"
)
# save computed indices
write.csv(ind$results, paste0(result_dir, "/11_sobol_indices.csv"), row.names = FALSE)

# Look at relevant parameters for ET
relevant_parameters <- c(# KSat
  "PTF_Ks_constant", "PTF_Ks_sand", "PTF_Ks_clay",
  # FieldCap
  "FieldCap_c1", "FieldCap_c2", # -0.6, 2.0 -> not optimized!
  # fRoots
  "rootFractionCoefficient_forest", "rootFractionCoefficient_impervious",
  "rootFractionCoefficient_pervious", 
  # ThetaS
  "PTF_lower66_5_constant", "PTF_lower66_5_clay", "PTF_lower66_5_Db",
  "PTF_higher66_5_constant", "PTF_higher66_5_clay", "PTF_higher66_5_Db",
  # L1_Max_Canopy_Intercept
  "canopyInterceptionFactor"
)
relevant_sobol_ind <- ind$results[ind$results$parameters %in% relevant_parameters, ]
write.csv(relevant_sobol_ind, paste0(result_dir, "/12_relevant_sobol_indices.csv"), row.names = FALSE)



















