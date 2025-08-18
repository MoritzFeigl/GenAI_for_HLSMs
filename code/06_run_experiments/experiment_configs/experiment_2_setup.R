# Experiment 2 setup
# optimization of 7 mHM TFs & numeric TF coefficients
# input data: standard mHM spatial inputs


# vaes ------------------------------------------------------------------------------------
# vae 1: base
# vae 2: base & base_plus
# vae 3: base & mhm
# vae 4: base & base_plus & mhm
# vae 5: lai
# Parameters: KSat, FieldCap, fRoots, ThetaS, L1_Max_Canopy_Intercept
tf_parameters <- 7
latent_dim <- 30
vae_list <- list("Ksat" = 2,
                 "FieldCap" = 4,
                 "ThetaS_1" = 2,
                 "ThetaS_2" = 2,
                 "fRoots_1" = 2,
                 "fRoots_2" = 2,
                 "L1_Max_Canopy_Intercept" = 5)


get_latent_bounds <- function(vae){
  bounds <- read.csv(paste0(data_dir, "08_trained_vaes/", "vae", vae, "_latent_bounds.csv"))
  bounds[, c("min", "max")]
}
GenAI_para_bounds <- lapply(vae_list, get_latent_bounds)
GenAI_para_bounds <- do.call(rbind, GenAI_para_bounds)

# scaling --------------------------------------------------------------------------------
scale_list <- list("FieldCap" = c(0.1, 1.0),
                   "ThetaS_1" = c(0.2, 1.0),
                   "ThetaS_2" = c(0.2, 1.0),
                   "fRoots_1" = c(0.1, 1.0),
                   "fRoots_2" = c(0.1, 1.0))

scale_factor <- "100.00"
# Numeric parameters ---------------------------------------------------------------------
# load mHM coefficients 
parameter_nml <- read_nml(paste0(data_dir, "05_mhm/master_mhm_config/mhm_parameter.nml"))
parameter_df <- do.call(rbind, parameter_nml$mhm_parameters)

# reduced parameter set defined by SA
red_standard_parameters <- readRDS(paste0(main_path, "results/04_SA/07_red_standard_parameters.rds"))
red_num_para_bounds <- readRDS(paste0(main_path, "results/04_SA/07_red_num_para_bounds.rds"))

# remove paras that are only relevant for KSat, FieldCap, ThetaS, fRoots, Canopy_Intercept
remove_old_tf_parameters <- c(# KSat
  "PTF_Ks_constant", "PTF_Ks_sand", "PTF_Ks_clay",
  # FieldCap
  "FieldCap_c1", "FieldCap_c2", # -0.6, 2.0 -> not optimized!
  # fRoots
  "rootFractionCoefficient_forest", "rootFractionCoefficient_pervious", 
  # ThetaS
  "PTF_lower66_5_constant", "PTF_lower66_5_clay", "PTF_lower66_5_Db",
  "PTF_higher66_5_constant", "PTF_higher66_5_clay", "PTF_higher66_5_Db",
  # L1_Max_Canopy_Intercept
  "canopyInterceptionFactor"
)
relevant_para_ids <- which(!(names(red_standard_parameters) %in% remove_old_tf_parameters))
GenAI_para_num_para <- red_standard_parameters[relevant_para_ids]
GenAI_para_num_para_bounds <- red_num_para_bounds[relevant_para_ids]

# get parameter indices and bounds
ind_of_num_paras <- which(rownames(parameter_df) %in% names(GenAI_para_num_para))
num_para_bounds <- parameter_df[ind_of_num_paras, c(1, 2)]
standard_parameters <- parameter_df[ind_of_num_paras, 3]

# SCE Parameters -------------------------------------------------------------------------
restart_after_n_iter <- 1
numIter <-  100
ncomplexes <- 2
xBounds <- data.frame(lower = c(GenAI_para_bounds$min-5, num_para_bounds[, 1]),
                      upper = c(GenAI_para_bounds$max+5, num_para_bounds[, 2]))
start_point <- c(rnorm(tf_parameters*latent_dim, 0, 1), standard_parameters)
NDIM <- length(start_point)


# Loss function mode ---------------------------------------------------------------------
# loss_modes:
# Q........ only discharge based loss
# QET...... discharge and ET based loss
# Q_TF..... discharge and TF based loss
# QET_TF... discharge, ET and TF based loss
loss_mode <- "Q_TF"