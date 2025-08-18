# Time mHM runs
# GenAI_para-ET-mHM project

# Run config
run_name <- "02_run_times"

# create run and result folders
run_folder <- paste0("/home/fs71468/mfeigl/GenAI_para_for_HLSMs/runs/", run_name)
result_dir <- paste0("/home/fs71468/mfeigl/GenAI_para_for_HLSMs/results/", run_name)
code_path <- "/home/fs71468/mfeigl/GenAI_para_for_HLSMs/code/"
data_path <- "/home/fs71468/mfeigl/GenAI_para_for_HLSMs/data/"

dir.create("/home/fs71468/mfeigl/GenAI_para_for_HLSMs/runs", showWarnings = FALSE)
dir.create(run_folder, showWarnings = FALSE)
dir.create("/home/fs71468/mfeigl/GenAI_para_for_HLSMs/results", showWarnings = FALSE)
dir.create(result_dir, showWarnings = FALSE)

setwd(run_folder)

# initial run directory creation
if(!dir.exists(paste0(run_folder, "/mhm"))){
  
  # mhm
  file.copy(paste0(data_path, "/05_mhm/mhm"), run_folder, 
            recursive=TRUE, overwrite = TRUE)
  
  # remove build folder
  try(unlink(paste0(run_folder, "/mhm/build"), recursive = TRUE))
  
  # copy and adapt first_run_preparation script
  file.copy(paste0(code_path, "02_utility_functions/04_first_run_preparation.sh"), 
            run_folder, overwrite = TRUE)
  run_preparation <- readLines(paste0(run_folder, "/04_first_run_preparation.sh"))
  run_preparation <- gsub("PROJECTPATH=/home/fs71468/mfeigl/GenAI_para_for_HLSMs/runs/",
                          paste0("PROJECTPATH=", run_folder), run_preparation)
  writeLines(run_preparation, paste0(run_folder, "/04_first_run_preparation.sh"))
  
  # copy and adapt run_mhm_preparation script
  file.copy(paste0(code_path, "02_utility_functions/05_run_mhm_preparation.sh"), 
            run_folder, overwrite = TRUE)
  run_preparation <- readLines(paste0(run_folder, "/05_run_mhm_preparation.sh"))
  run_preparation <- gsub("PROJECTPATH=/home/fs71468/mfeigl/GenAI_para_for_HLSMs/runs/",
                          paste0("PROJECTPATH=", run_folder), run_preparation)
  writeLines(run_preparation, paste0(run_folder, "/05_run_mhm_preparation.sh"))
  
  # create output directories
  dir.create(paste0(run_folder, "/Training"), showWarnings = FALSE)
  
  # adapt paths in mhm.nml for each basin and change time period
  basins <- list.files(paste0(data_path, "06_training_validation_data/Training/"))
  
  for(basin in basins){
    # create basin specific mhm config and output folder
    dir.create(paste0(run_folder, "/Training/", basin), showWarnings = FALSE)
    basin_dir <- paste0(run_folder, "/Training/", basin)
    file.copy(paste0(data_path, "06_training_validation_data/Training/", basin, "/config/"),
      basin_dir, recursive=TRUE, overwrite = TRUE)
    
    
    # mhm.nml
    file.copy(paste0(data_path, "05_mhm/master_mhm_config/", "/mhm.nml"),
              paste0(basin_dir,"/config/mhm.nml"), overwrite = TRUE)
    mhm_nml <- readLines(paste0(basin_dir,"/config/mhm.nml"))
    mhm_nml <- sapply(mhm_nml, function(x) gsub("<SPLIT>", "Training", x))
    mhm_nml <- sapply(mhm_nml, function(x) gsub("<BASIN>", basin, x))
    writeLines(mhm_nml, paste0(basin_dir,"/config/mhm.nml"))
    
    # mpr.nml
    file.copy(paste0(data_path, "05_mhm/master_mhm_config/", "mpr.nml"),
              paste0(basin_dir,"/config/mpr.nml"), overwrite = TRUE)
    mpr_nml <- readLines(paste0(basin_dir,"/config/mpr.nml"))
    mpr_nml <- sapply(mpr_nml, function(x) gsub("<SPLIT>", "Training", x))
    mpr_nml <- sapply(mpr_nml, function(x) gsub("<BASIN>", basin, x))
    writeLines(mpr_nml, paste0(basin_dir,"/config/mpr.nml"))
    
    file.copy(paste0(data_path, "/05_mhm/master_mhm_config/run_single_basin.sh"),
              paste0(basin_dir, "/config"), recursive=TRUE, overwrite = TRUE)
    
    file.copy(paste0(data_path, "06_training_validation_data/Training/", basin, "/output"),
              basin_dir, recursive=TRUE, overwrite = TRUE)
  }
  
  # initial compilation
  system("bash 04_first_run_preparation.sh & wait")
  
}

basins <- list.files(paste0(run_folder, "/Training"))
run_times <- data.frame(basin = basins, start = NA, end = NA, time = NA)
for(i in seq_along(basins)){
  cat("Run basin", basins[i], paste0(i, "/", length(basins)), "\n")
  start <- Sys.time()
  setwd(paste0("Training/", basins[i], "/config"))
  system("bash run_single_basin.sh")
  end <- Sys.time()
  run_times$start[i] <- as.numeric(start)
  run_times$end[i] <- as.numeric(end)
  run_times$time[i] <- as.numeric(end) - as.numeric(start)
  setwd(run_folder)
  write.csv(run_times, paste0(result_dir, "/training_basin_run_times.csv"))
}

