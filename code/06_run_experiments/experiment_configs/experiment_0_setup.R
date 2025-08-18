# Experiment 0 setup
# optimization of standard mHM TF coefficients

# VAEs -----------------------------------------------------------------------------------
# None

# Numeric parameters ---------------------------------------------------------------------
# load mHM coefficients 
parameter_nml <- read_nml(paste0(data_dir, "05_mhm/master_mhm_config/mhm_parameter.nml"))
parameter_df <- do.call(rbind, parameter_nml$mhm_parameters)

# reduced parameter set defined by SA
red_standard_parameters <- readRDS(paste0(main_path, "results/04_SA/07_red_standard_parameters.rds"))
red_num_para_bounds <- readRDS(paste0(main_path, "results/04_SA/07_red_num_para_bounds.rds"))

# get parameter indices and bounds
ind_of_num_paras <- which(rownames(parameter_df) %in% names(red_standard_parameters))
num_para_bounds <- parameter_df[ind_of_num_paras, c(1, 2)]
standard_parameters <- parameter_df[ind_of_num_paras, 3]
# SCE Parameters -------------------------------------------------------------------------
restart_after_n_iter <- 2
numIter <-  100
ncomplexes <- 2
xBounds <- data.frame(lower = num_para_bounds[, 1],
                      upper = num_para_bounds[, 2])
start_point <- standard_parameters
NDIM <- length(start_point)

# Loss function mode ---------------------------------------------------------------------
# loss_modes:
# Q........ only discharge based loss
# QET...... discharge and ET based loss
# Q_TF..... discharge and TF based loss
# QET_TF... discharge, ET and TF based loss
loss_mode <- "Q"












