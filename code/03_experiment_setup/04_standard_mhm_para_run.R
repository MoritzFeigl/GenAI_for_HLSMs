# Time mHM runs
# GenAI_para-ET-mHM project

# Run config
run_name <- "03_standard_mhm_para_run"
loss_mode = "QET"

# create run and result folders
run_dir <- paste0("/gpfs/data/fs71468/GenAI_para_runs/runs/", run_name)
result_dir <- paste0("/gpfs/data/fs71468/GenAI_para_runs/results/", run_name)
code_dir <- "/gpfs/data/fs71468/GenAI_para_runs/code/"
data_dir <- "/gpfs/data/fs71468/GenAI_para_runs/data/"
dir.create("/gpfs/data/fs71468/GenAI_para_runs/runs", showWarnings = FALSE)
dir.create(run_dir, showWarnings = FALSE)
dir.create("/gpfs/data/fs71468/GenAI_para_runs/results", showWarnings = FALSE)
dir.create(result_dir, showWarnings = FALSE)

# set wd and source utils
setwd(run_dir)
source(paste0(code_dir, "02_utility_functions/01_GenAI_para_utility_functions.R"))

# create run environment in run_dir
create_run_environment(run_dir, data_dir, code_dir, n_parallel_setups=1)

# Change output to only aET monthly
state_and_fluxes[["aET"]] <- TRUE
change_mhm_outputs(basins_dir="Training_1", state_and_fluxes, time_step = "monthly")

# get numeric parameters
parameter_nml <- nml::read_nml(paste0(data_dir, "05_mhm/master_mhm_config/mhm_parameter.nml"))
parameter_df <- do.call(rbind, parameter_nml$mhm_parameters)
ind_of_num_paras <- which(parameter_df[, 4] == 1)
num_para_bounds <- parameter_df[ind_of_num_paras, c(1, 2)]
standard_parameters <- parameter_df[ind_of_num_paras, 3]

# Standard parameter
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
update_tfs_in_mpr_nml(run_dir=run_dir, setup = 1, 
                      parameter_list=parameter_list, scaling = FALSE, save_para_as_ncdf = TRUE)
update_numeric_parameters(numeric_parameters = standard_parameters, 
                          run_dir = run_dir, setup=1, parameter_nml=parameter_nml, 
                          ind_of_num_paras=ind_of_num_paras)

# Runs mHM
start_eval <- Sys.time()
run_mhm(run_dir)
end_eval <- Sys.time()
cat("Evaluation took:\n")
print(end_eval - start_eval)

# get ET and Q results and compute quality criteria
basin_qc <- basins_quality_criteria(run_dir, setup, data_dir)
loss_list <- compute_loss(basin_qc, loss_mode, point_tf)
result <- cbind(basin_qc, matrix(loss_list[[2]], ncol=1))
names(result)[ncol(result)] <- paste0(loss_mode, "-loss")
write.csv(result, paste0(result_dir, "/basin_results.csv"))


