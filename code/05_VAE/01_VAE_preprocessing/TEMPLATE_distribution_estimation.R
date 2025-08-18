# Estimate TF distributions for GenAI_para
# GenAI_para-ET-mHM project

main_path <- "/gpfs/data/fs71468/GenAI_para_runs/"
code_dir <- paste0(main_path, "code")
code_files <- list.files(paste0(code_dir, "/02_utility_functions/09_CFG_implementation"), full.names = TRUE)
for(file in code_files) source(file)
setwd(paste0(main_path, "data/07_VAE_data/03_prepared_TFs/prepared_functions"))

# Start script
variable_df_sampled <- feather::read_feather(paste0(main_path, "data/07_VAE_data/01_sp_data/spatial_predictors_sampled.feather"))
variable_df_sampled <- variable_df_sampled[, -which(names(variable_df_sampled) == "basin")]
# Use only till variables for distribution estimation
variable_df_sampled$KSat_notill <- NULL
variable_df_sampled$vGenu_n_notill <- NULL
variable_df_sampled$ThetaS_notill <- NULL
names(variable_df_sampled) <- gsub("_till", "", names(variable_df_sampled))

# add aspect sin/cos transformation
variable_df_sampled$aspect_sin <- sin(variable_df_sampled$aspect*(2*pi/360))
variable_df_sampled$aspect_cos <- cos(variable_df_sampled$aspect*(2*pi/360))
variable_df_sampled$aspect <- NULL

# Compute function distributions
functions_for_dist <- "dummy_file"

# # scaling_bounds
# scaling_bounds <- list("slope" = c(0, 90),
#                        "aspect_sin" = c(-1, 1),
#                        "aspect_cos" = c(-1, 1),
#                        "bd" = c(1.1, 2.3), #max(variable_df$bd)), # is actually 2.2769
#                        "sand" = c(0, 100),
#                        "clay" = c(0, 100),
#                        "dem" = c(0, 2000),
#                        "KSat" = c(1.1, 1000),
#                        "vGenu_n" = c(1, 2),
#                        "ThetaS" = c(0.01, 0.55),
#                        "mat" = c(0, 20),
#                        "mat_range" = c(15, 25),
#                        "map" = c(200, 3000),
#                        "lai" = c(0, 7)
# )
# calc dist for all given files
functions_para <- fst::read_fst(functions_for_dist)
dist_name <- gsub("simplified", "distribution", functions_for_dist)
distribution_sampler(functions = functions_para,
                     variable_df = variable_df_sampled,
                     scaling_bounds = FALSE,
                     file_name = dist_name,
                     no_cores = 120)
