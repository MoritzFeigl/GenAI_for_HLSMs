# GenAI_para-ET-mHM project
# analyze validation

library(magrittr)
library(ggplot2)
library(wesanderson)
run_name <- "Validation"
main_path <- "/gpfs/data/fs71468/GenAI_para_runs/"
training_results_file <- paste0(main_path, "results/training_results_2024-01-16.csv")

# Create validation environment ----------------------------------------------------------
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
basins <- list.files(paste0(data_dir, "06_training_validation_data/Validation"))
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
  run_preparation <- gsub("/Training", paste0("/Validation_", setup), run_preparation, fixed = TRUE)
  run_preparation <- gsub("6335125", basins[1], run_preparation)
  run_preparation <- gsub("<DATA_DIR>", data_dir, run_preparation, fixed = TRUE)
  writeLines(run_preparation, paste0(run_dir, "/04_first_run_preparation_", setup, ".sh"))
  
  # run preparation (mhm compilation) script
  run_preparation <- readLines(paste0(code_dir, "02_utility_functions/05_run_mhm_preparation.sh"))
  run_preparation <- gsub("PROJECTPATH=/home/fs71468/mfeigl/GenAI_para_for_HLSMs/runs/",
                          paste0("PROJECTPATH=", run_dir), run_preparation, fixed = TRUE)
  run_preparation <- gsub("/mhm/", paste0("/mhm_", setup, "/"), run_preparation, fixed = TRUE)
  run_preparation <- gsub("/Training", paste0("/Validation_", setup), run_preparation, fixed = TRUE)
  run_preparation <- gsub("6335125", basins[1], run_preparation, fixed = TRUE)
  run_preparation <- gsub("<DATA_DIR>", data_dir, run_preparation, fixed = TRUE)
  writeLines(run_preparation, paste0(run_dir, "/05_run_mhm_preparation_", setup, ".sh"))
  
  # create basin mhm dirs
  dir.create(paste0(run_dir, "/Validation_", setup), showWarnings = FALSE)
  
  # adapt paths in mhm.nml for each basin and change time period
  for(basin in basins){
    # create basin specific mhm config and output folder
    basin_dir <- paste0(run_dir, "/Validation_", setup, "/", basin)
    dir.create(basin_dir, showWarnings = FALSE)
    file.copy(paste0(data_dir, "06_training_validation_data/Validation/", basin, "/config/"),
              basin_dir, recursive=TRUE, overwrite = TRUE)
    
    # mhm.nml
    mhm_nml <- readLines(paste0(basin_dir,"/config/mhm.nml"))
    mhm_nml <- sapply(mhm_nml, function(x) gsub("<DATA_DIR>", data_dir, x))
    mhm_nml <- sapply(mhm_nml, function(x) gsub("<SPLIT>", "Validation", x))
    mhm_nml <- sapply(mhm_nml, function(x) gsub("<BASIN>", basin, x))
    mhm_nml <- sapply(mhm_nml, function(x) gsub("/Training/", "/Validation/", x, fixed = TRUE))
    writeLines(mhm_nml, paste0(basin_dir,"/config/mhm.nml"))
    
    # mpr.nml
    file.copy(paste0(data_dir, "05_mhm/master_mhm_config/", "mpr.nml"),
              paste0(basin_dir,"/config/mpr.nml"), overwrite = TRUE)
    mpr_nml <- readLines(paste0(basin_dir,"/config/mpr.nml"))
    mpr_nml <- sapply(mpr_nml, function(x) gsub("<DATA_DIR>", data_dir, x,))
    mpr_nml <- sapply(mpr_nml, function(x) gsub("<SPLIT>", "Validation", x))
    mpr_nml <- sapply(mpr_nml, function(x) gsub("/Training/", "/Validation/", x, fixed = TRUE))
    mpr_nml <- sapply(mpr_nml, function(x) gsub("<BASIN>", basin, x))
    writeLines(mpr_nml, paste0(basin_dir,"/config/mpr.nml"))
    
    # create basin output directory
    file.copy(paste0(data_dir, "06_training_validation_data/Validation/", basin, "/output"),
              basin_dir, recursive=TRUE, overwrite = TRUE)
    
  }
  
  # Create run scripts
  dir.create(paste0("Validation_run_scripts_", setup), showWarnings = FALSE)
  for(basin in list.files(paste0(run_dir, "/Validation_", setup))){
    run_script_name <- paste0(run_dir, "/Validation_run_scripts_", setup, 
                              "/basin_", basin, ".sh")
    fileConn <- file(run_script_name)
    writeLines(
      c("#!/bin/bash",
        paste0("cd ", run_dir, "/Validation_", setup, "/", basin, "/config"),
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

# 
# # find GenAI_para runs that produced valid GenAI_para functions
# training_results <- read.csv(training_results_file)
# 
# validity <- data.frame("experiment"=rep(1:5, each=5), "run"=1:5, "valid"=NA, "loss"=NA)
# for(experiment in 1:4){
#   for(run in 1:5){
#     try({
#       exp_results <- unique(training_results[training_results$experiment == experiment & 
#                                                training_results$run == run, c("run", "mean_loss", "time")])
#       # Find real best result if VAR_NUM_ERROR
#       variables <- c("dem", "aspect_sin", "aspect_cos", "slope", "bd", "sand", "clay", "lai", "map", "mat", "mat_range", "ThetaS", "KSat", "vGenu_n", "lai")
#       tf_tracker <- rbind(read.csv(paste0("/gpfs/data/fs71468/GenAI_para_runs/results/", "exp", experiment, "-run", run, "/TFs_tracker_1.csv")),
#                           read.csv(paste0("/gpfs/data/fs71468/GenAI_para_runs/results/", "exp", experiment, "-run", run, "/TFs_tracker_2.csv")))
#       
#       best_id <- which(tf_tracker$time == exp_results$time)
#       print(best_id)
#       
#       if(sum(is.na(as.numeric(sapply(tf_tracker[best_id, 4:10], tf_evaluation, variables=variables)))) == 0){
#         validity[validity$experiment == experiment & validity$run == run, "valid"] <- TRUE
#         validity[validity$experiment == experiment & validity$run == run, "loss"] <- exp_results$mean_loss
#       } else {
#         validity[validity$experiment == experiment & validity$run == run, "valid"] <- FALSE
#       }
#     })
#   }
# }

# Compute validation results -------------------------------------------------------------

for(experiment in c(0, 1, 2)){
  
  training_results <- read.csv(training_results_file)
  exp_results <- unique(training_results[training_results$experiment == experiment, 
                                         c("run", "mean_loss", "time")])
  run <- exp_results[exp_results$mean_loss == max(exp_results$mean_loss), "run"]
  
  # # Find real best result if VAR_NUM_ERROR
  # variables <- c("dem", "aspect_sin", "aspect_cos", "slope", "bd", "sand", "clay", "lai", "map", "mat", "mat_range", "ThetaS", "KSat", "vGenu_n", "lai")
  # tf_tracker <- read.csv(paste0("/gpfs/data/fs71468/GenAI_para_runs/results/", "exp", experiment, "-run", run, "/TFs_tracker_1.csv"))
  # 
  # best_id <- which(tf_tracker$time == exp_results$time[run])
  # tf_tracker[best_id, ]
  # checked_tfs <- NULL
  # 
  # for(tf in c(names(tf_tracker)[4:10])){
  #   print(tf)
  #   checked_tfs <- cbind(checked_tfs, as.numeric(sapply(tf_tracker[, tf], tf_evaluation, variables=variables)))
  # }
  # 
  # for(i in best_id:1){
  #   cat(i, "\r")
  #   if(sum(is.na(checked_tfs[i, ])) == 0) break
  # }
  
  # define experiment configuration
  source(paste0(code_dir, 
                "06_run_experiments/experiment_configs/experiment_", 
                experiment, "_setup.R"))
  
  # define mhm output --> always aET for result analysis
  state_and_fluxes[["aET"]] <- TRUE
  for(setup in 1:n_parallel_setups){
    change_mhm_outputs(basins_dir = paste0(run_dir, "/Validation_", setup), 
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
  
  # custom
  # KSat <- "-(dem*sqrt(mat-sand+66.76*slope-abs(map)))-sand"
  # FieldCap <- "( 8.43/abs((map)-439.44+vGenu_n+mat) ) / 100.00"
  # fRoots_1 <- "( abs(dem-38296/646.07^-26.92-bd7)/mat_range ) / 100.00"
  # fRoots_2 <- "( tan(27771.52)/(mat+49.87)-slope/29.12*sand ) / 100.00"
  # ThetaS_1 <- "( -(abs(sand)-mat_range/sin(sand93)-95.53/dem) ) / 100.00"
  # ThetaS_2 <- "( 356.43-14.98/sand ) / 100.00"
  # Canopy_Intercept <- "0.76/2.55/lai+4"
  
  parameter_list <- list("Ksat" = KSat,
                         "FieldCap" = FieldCap,
                         "fRoots_1" = fRoots_1, "fRoots_2" = fRoots_2,
                         "ThetaS_1" = ThetaS_1, "ThetaS_2" = ThetaS_2,
                         "L1_Max_Canopy_Intercept" = Canopy_Intercept)
  
  # if all TFs are valid run mHM and compute losses
  update_tfs_in_mpr_nml(run_dir, setup, parameter_list, scaling = FALSE, 
                        save_para_as_ncdf = TRUE, Training=FALSE)
  
  # select numeric mhm para based on best results
  name_fix_df <- data.frame(matrix(NA, nrow=1, ncol=length(names(standard_parameters))))
  names(name_fix_df) <- names(standard_parameters)
  relevant_parameters <- names(data.frame(name_fix_df))
  best_num_para <- best_parameter[, relevant_parameters]
  
  # numeric para update
  update_numeric_parameters(best_num_para, run_dir, setup, 
                            parameter_nml, ind_of_num_paras, Training=FALSE)
  
  for(basin in basins){
    basin_dir <- paste0(run_dir, "/Validation_", setup, "/", basin)
    # mpr.nml
    mpr_nml <- readLines(paste0(basin_dir,"/config/mpr.nml"))
    mpr_nml <- sapply(mpr_nml, function(x) gsub("/Training/", "/Validation/", x, fixed = TRUE))
    writeLines(mpr_nml, paste0(basin_dir,"/config/mpr.nml"))
    mhm_nml <- readLines(paste0(basin_dir,"/config/mhm.nml"))
    mhm_nml <- sapply(mhm_nml, function(x) gsub("/Training/", "/Validation/", x, fixed = TRUE))
    writeLines(mhm_nml, paste0(basin_dir,"/config/mhm.nml"))
  }
  
  # Run mHM
  run_mhm(run_dir, setup, mhm_preparation = TRUE, Training=FALSE)
  
  # # run single basins
  # basin1 <- "6340300"
  # basin2 <- "9304066"
  # system(paste0("bash 05_run_mhm_preparation_", setup, ".sh & wait"))
  # basin = "9316288"
  # basin_dir <- paste0(run_dir, "/Validation_", setup, "/", basin)
  # 
  # setwd(paste0(basin_dir, "/config"))
  # system("export LD_LIBRARY_PATH=$LIBRARY_PATH
  # ./mhm")
  # system(paste0("export LD_LIBRARY_PATH=$LIBRARY_PATH
  # cd Validation_run_scripts_", setup, "
  #               bash basin_", basin1, ".sh &
  #               bash basin_", basin2, ".sh & wait"))
  
  # compute quality criteria and loss
  basin_qc <- basins_quality_criteria(run_dir, setup, data_dir, Training=FALSE, ET=TRUE)
  
  # save results, parameter fields and fluxes
  exp_result_dir <- paste0(main_path, "results/", run_name, "/exp", experiment, "-run", run, "/")
  dir.create(exp_result_dir)
  
  # 1 save basin_qc and loss
  write.csv(basin_qc, 
            paste0(exp_result_dir, "basins_quality_criteria_validation.csv"), 
            row.names = FALSE)
  
  # # 2. save parameter files
  # for(basin in basins){
  #   system(paste0("cp -f ", 
  #                 paste0(run_dir, "/Validation_1/", basin, "/output/mHM_parameters.nc "),
  #                 paste0(exp_result_dir, basin, "_parameters.nc")))
  # }
  # 
  # # 3. save fluxes
  # for(basin in basins){
  #   system(paste0("cp -f ", 
  #                 paste0(run_dir, "/Validation_1/", basin, "/output/mHM_Fluxes_States.nc "),
  #                 paste0(exp_result_dir, basin, "_Fluxes_States.nc")))
  # }
  # 
  # 4. save discharge
  for(basin in basins){
    system(paste0("cp -f ", 
                  paste0(run_dir, "/Validation_1/", basin, "/output/daily_discharge.out "),
                  paste0(exp_result_dir, basin, "_daily_discharge.out")))
  }
  
}
# Plot QC Boxplots -----------------------------------------------------------------------

all_qc <- NULL
for(exp_run in list.files(result_dir)){
  qc <- read.csv(paste0(result_dir, "/", exp_run, "/basins_quality_criteria_validation.csv"))
  all_qc <- rbind(all_qc, cbind("experiment"=as.integer(substr(exp_run, 4, 4)), 
                                "run"=as.integer(substr(exp_run, 9, 9)), qc))
}

plot_data <- reshape2::melt(all_qc, id.vars=c("experiment", "run", "Basin"))

# function to get median values in boxplots
get_box_stats <- function(y, upper_limit=1) {
  return(data.frame(
    y = upper_limit,
    label = paste(
      "Median =", round(median(y), 2)
    )
  ))
}

plot_data$tf_type <- "GenAI_para"
plot_data$tf_type[plot_data$experiment == 0] <- "mHM"
plot_data$tf_type <- factor(plot_data$tf_type, levels = c("GenAI_para", "mHM"))

plot_data$experiment <- factor(plot_data$experiment, levels = c(0, 1, 2),
                               labels = c("mHM", "GenAI_para", "GenAI_para + SP"))
plot_data$variable <- factor(plot_data$variable, levels=c("NSE", "lNSE", "KGE", "SPAEF"))

ggplot(plot_data, 
       aes(factor(experiment), value, fill=tf_type)) + 
  geom_boxplot() + 
  geom_jitter(width = 0.2, alpha=0.4) + 
  facet_wrap(~variable, scales = "free")  + 
  coord_cartesian(ylim = c(0,1)) +
  scale_fill_manual("TF type", values=wes_palette(name="Darjeeling1", n = 4)[c(2,4)]) +
  labs(x="Experiment", fill="TF type", y="") +
  stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0, size=3)
ggsave(paste0(main_path, "/analysis_results/validation_results/best_validation_runs.png"),
       width = 18, height = 8, units = "in")

experiment_subset <- c("mHM", "GenAI_para", "GenAI_para + SP")
ggplot(plot_data[plot_data$variable != "SPAEF" &
                   plot_data$experiment %in% experiment_subset, ], 
       aes(factor(experiment), value, fill=tf_type)) + 
  geom_boxplot() + 
  geom_jitter(width = 0.2, alpha=0.4) + 
  facet_wrap(~variable, scales = "free")  + 
  coord_cartesian(ylim = c(0,1)) +
  scale_fill_manual("TF type", values=wes_palette(name="Darjeeling1", n = 4)[c(2,4)]) +
  labs(x="Experiment", fill="TF type", y="") +
  stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0, size=3)
ggsave(paste0(main_path, "/analysis_results/validation_results/best_validation_runs_NoSpaef.png"),
       width = 18, height = 8, units = "in")


write.csv(plot_data, paste0(main_path, "/analysis_results/validation_results/validation_plot_data.csv"))


# Combine with basin info
# main_path <- "C:/Users/morit/Dropbox/Projekte/GenAI_para_for_HLSMs"
# plot_data <- read.csv(paste0(main_path, "/analysis_results/validation_results/validation_plot_data.csv"))
# 
# basin_info <- readxl::read_excel(paste0(main_path, "/data/01_study_basins_info/basins_master_file.xlsx"), 1)
# 
# plot_data_info <- merge(plot_data, basin_info, by.x="Basin", by.y="Stat_ID")
# 
# for(experiment in unique(plot_data$experiment)){
#   for(var in c("NSE", "lNSE", "KGE")){
#   sub_data <- plot_data_info[plot_data_info$experiment == experiment & 
#                                plot_data_info$variable == var, ]
#   names(sub_data)[30] <- "Qobs"
#   sum_res <- summary(lm(formula("value~catArea + dem_mean + Qobs"), sub_data))
#   if(any(sum_res$coefficients[-1, 4] < 0.05)){
#     cat(experiment, "-", var, ":")
#     sum_res$coefficients[, 4] <- round(sum_res$coefficients[, 4], 3)
#     print(sum_res$coefficients[, c(1, 4)])
#     cat("\n\n")
#   }
#   }
# }

