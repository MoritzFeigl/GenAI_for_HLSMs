# GenAI_para_for_HLSMs utility functions

if(!("ncdf4" %in% installed.packages())){
  chooseCRANmirror(ind = 6)
  install.packages("ncdf4")
}
suppressWarnings(library(ncdf4))

if(!("feather" %in% installed.packages())){
  chooseCRANmirror(ind = 6)
  install.packages("feather")
}

if(!("magrittr" %in% installed.packages())){
  chooseCRANmirror(ind = 6)
}
suppressWarnings(library(magrittr))


# TF utilities ---------------------------------------------------------------------------
.function_splitter <- function(point_tf){
  function_splitted <- unlist(strsplit(point_tf, c("[/^()*+-]")))
  function_splitted <- gsub(" ", "", function_splitted)
  function_splitted <- function_splitted[function_splitted != ""]
  return(function_splitted)
}

# function that adds correct scaling numerical values to tf
prepare_tf_for_mhm <- function(tf, scaling, scaling_bounds = NULL, variables){
  
  if(scaling){
    for(i in 1:length(scaling_bounds)){
      scaling_bounds[[i]] <- scaling_bounds[[i]] %>% format(nsmall = 1) %>% gsub(" ", "", .)
    }
  }
  # prepare
  tf <- gsub(" ", "", tf)
  tf <- gsub("^", "**", tf, fixed = TRUE)
  # split tf
  tf_splitted <- tf %>%
    strsplit(., c("[/^()*+-]")) %>%
    unlist()
  var <- character()
  for(i in seq_along(tf_splitted)){
    if(tf_splitted[i] %in% variables)  var <- c(var, tf_splitted[i])
  }
  var <- unique(var)
  
  if(scaling){
    # change variables to scaled variables in tf
    for(variable in var){
      bounds <- scaling_bounds[[variable]]
      tf <- gsub(pattern = variable,
                 replacement = paste0("((", variable, "-", bounds[1], ")/(",
                                      bounds[2], "-", bounds[1], "))"),
                 tf)
    }
  }
  # get numerical values
  tf_splitted <- tf %>%
    strsplit(., c("[/^()*+-]")) %>%
    unlist()
  #tf_splitted <- tf_splitted[tf_splitted != ""]
  numerics <- character()
  suppressWarnings(
    for(i in seq_along(tf_splitted)){
      id_num <- ifelse(is.na(as.numeric(tf_splitted[i])), NA, tf_splitted[i])
      numerics <- c(numerics, id_num)
    }
  )
  
  tf_splitted_full <- tf %>%
    strsplit(., "(?<=[/^()*+-])", perl=TRUE) %>%
    unlist()
  
  
  # make integers to float
  ints <- grep(".", numerics, invert = TRUE, fixed = TRUE)
  
  if(length(ints) != 0){
    for(int in ints) {
      if(!is.na(numerics[int])){
        tf_splitted_full[int] <- gsub(numerics[int], paste0(numerics[int], ".00"), 
                                      tf_splitted_full[int])
        numerics[int] <- paste0(numerics[int], ".00")
      }
    }
  }
  
  # reconstruct tf
  tf <- paste0(tf_splitted_full, collapse = "")
  
  # remove numeric NAs
  numerics <- numerics[!is.na(numerics)]
  
  # put space into function
  tf <- tf %>%
    gsub("+", " + ", ., fixed = TRUE) %>%
    gsub("-", " - ", ., fixed = TRUE) %>%
    gsub("*", " * ", ., fixed = TRUE) %>%
    gsub("*  *", "**", ., fixed = TRUE) %>%
    gsub("/", " / ", ., fixed = TRUE)
  
  if(substr(tf, 1, 1) == " ") tf <- substring(tf, 2)
  return(list("function" = tf, "numerics" = numerics))
}


# Optimization utilities -----------------------------------------------------------------
state_and_fluxes <- list("L1_inter"=FALSE, 
                         "L1_snowpack"=FALSE, 
                         "L1_soilMoist"=FALSE, 
                         "vol_soilMoist_layer"=FALSE,
                         "vol_soilMoist_average"=FALSE,
                         "L1_sealSTW"=FALSE,
                         "L1_unsatSTW"=FALSE,
                         "L1_satSTW"=FALSE,
                         "PET"=FALSE,
                         "aET"=TRUE,
                         "L1_total_runoff"=FALSE,
                         "L1_runoffSeal"=FALSE,
                         "L1_fastRunoff"=FALSE,
                         "L1_slowRunoff"=FALSE,
                         "L1_baseflow"=FALSE,
                         "L1_percol"=FALSE,
                         "L1_infilSoil"=FALSE,
                         "ground_abedo"=FALSE,
                         "soil_evap"=FALSE,
                         "L1_preEffect"=FALSE)

change_mhm_outputs <- function(basins_dir, state_and_fluxes, time_step = "monthly"){
  # get basin folders
  basins <- list.files(basins_dir)
  # set types of mhm outputs
  mhm_out <- readLines(paste0(basins_dir, "/", basins[1], "/config/mhm_outputs.nml"))
  for(i in 1:length(state_and_fluxes)){
    id <- grep(paste0("(", i, ")"), mhm_out, fixed = TRUE)
    mhm_out[id] <- gsub(as.character(!state_and_fluxes[[i]]), 
                        as.character(state_and_fluxes[[i]]), 
                        mhm_out[id])
  }
  # set time_step
  time_setting <- switch(time_step,
                         at_end = 0,
                         daily = -1,
                         monthly = -2,
                         yearly = -3)
  id <- grep("timeStep_model_outputs", mhm_out, fixed = TRUE)
  mhm_out[id] <- paste0("timeStep_model_outputs = ", time_setting)
  for(basin in basins){
    writeLines(mhm_out, con = paste0(basins_dir, "/", basin, "/config/mhm_outputs.nml"))
  }
}

create_run_environment <- function(run_dir, data_dir, code_dir, n_parallel_setups = 4){
  
  # Check if environment exists
  if(dir.exists(paste0(run_dir, "/training_run_scripts_", n_parallel_setups))){
    return(cat("Environment exists already.\n"))
  }
  
  for(setup in 1:n_parallel_setups){
    
    # basin list
    basins <- list.files(paste0(data_dir, "06_training_validation_data/Training"))
    
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
    run_preparation <- gsub("/Training", paste0("/Training_", setup), run_preparation, fixed = TRUE)
    run_preparation <- gsub("6335125", basins[1], run_preparation)
    run_preparation <- gsub("<DATA_DIR>", data_dir, run_preparation, fixed = TRUE)
    writeLines(run_preparation, paste0(run_dir, "/04_first_run_preparation_", setup, ".sh"))
    
    # run preparation (mhm compilation) script
    run_preparation <- readLines(paste0(code_dir, "02_utility_functions/05_run_mhm_preparation.sh"))
    run_preparation <- gsub("PROJECTPATH=/home/fs71468/mfeigl/GenAI_para_for_HLSMs/runs/",
                            paste0("PROJECTPATH=", run_dir), run_preparation, fixed = TRUE)
    run_preparation <- gsub("/mhm/", paste0("/mhm_", setup, "/"), run_preparation, fixed = TRUE)
    run_preparation <- gsub("/Training", paste0("/Training_", setup), run_preparation, fixed = TRUE)
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
      mpr_nml <- sapply(mpr_nml, function(x) gsub("<DATA_DIR>", data_dir, x))
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
  }
  # # initial compilation
  comp_cmd <- paste0(paste0(
    "bash 04_first_run_preparation_",
    1:n_parallel_setups, ".sh & ", collapse = ""), "wait")
  
  system(comp_cmd)
}

# create pre-defined training batchelss for multi node computation
create_batches_for_multi_node_mhm <- function(run_dir, nr_nodes, run_times_path){
  # check if batches already exists
  if(dir.exists(paste0(run_dir, "/training_run_scripts"))){
    if(length(list.files(paste0(run_dir, "/training_run_scripts"))) == nr_nodes){
      return(cat("The batch run files are already prepared."))
    }
  }
  # get run times of all traiing basins
  run_times <- read.csv(run_times_path)
  # define batches for nr_nodes
  assign.job <- function(machines, job) {
    which.machines <- which.min(lapply(machines, sum))
    machines[[which.machines]] <- c(machines[[which.machines]], job)
    machines
  }
  allocate <- function(num.machines, job.times) {
    machines <- lapply(1:num.machines, function(...) c())
    Reduce(assign.job,
           sort(job.times, decreasing=TRUE),
           machines)
  }
  groups <- allocate(nr_nodes, run_times$time)
  batches <- NULL
  for (batch in 1:nr_nodes){
    batches <- rbind(batches, data.frame(batch = batch, time = groups[[batch]]))
  }
  batch_df <- merge(run_times, batches, by = "time")
  training_batches <- data.frame(basin = batch_df$basin, 
                                 time = batch_df$time, 
                                 batch = batch_df$batch)
  write.csv(training_batches, 
            paste0(run_dir, "/training_basins_batches.csv"),
            row.names = FALSE)
  
  # Create bash files for starting nodes
  unlink("training_run_scripts", recursive=TRUE)
  dir.create("training_run_scripts", showWarnings = FALSE)
  
  # create run script in train_run_scripts dir
  for(basin in list.files(paste0(run_dir, "/Training"))){
    basin_batch <- training_batches[training_batches$basin == basin, "batch"]
    dir.create(paste0(run_dir, "/training_run_scripts/batch_", basin_batch), showWarnings = FALSE)
    
    run_script_name <- paste0(run_dir, "/training_run_scripts/batch_",
                              basin_batch, "/run_" ,basin, ".sh")
    fileConn <- file(run_script_name)
    writeLines(
      c("#!/bin/bash",
        paste0("cd ", run_dir, "/Training/", basin, "/config"),
        "./mhm > output.txt"),
      fileConn)
    close(fileConn)
  }
  # create batch script in train_run_scripts dir
  for(batch in unique(training_batches$batch)){
    fileConn <- file(paste0(run_dir, "/training_run_scripts/run_batch_", batch, ".sh"))
    basin_sh <- list.files(paste0(run_dir, "/training_run_scripts/batch_", batch))
    
    
    writeLines(
      c("#!/bin/sh",
        "source ~/env/spackenv",
        paste0(
          paste0("bash ./batch_", batch, "/", basin_sh[-length(basin_sh)], collapse = " & "),
          " & process_id=$! bash ./batch_", batch, "/", basin_sh[length(basin_sh)], " & wait $process_id"
        ),
        'date +"%Y-%m-%d %T"'),
      fileConn)
    close(fileConn)
  }
  return(cat("Finished preparing all batch run files."))
}

# functions that makes all console output from system call available
robust.system <- function (cmd) {
  stderrFile = tempfile(pattern="R_robust.system_stderr", fileext=as.character(Sys.getpid()))
  stdoutFile = tempfile(pattern="R_robust.system_stdout", fileext=as.character(Sys.getpid()))
  
  retval = list()
  retval$exitStatus = system(paste0(cmd, " 2> ", shQuote(stderrFile), " > ", shQuote(stdoutFile)))
  retval$stdout = readLines(stdoutFile)
  retval$stderr = readLines(stderrFile)
  
  unlink(c(stdoutFile, stderrFile))
  return(retval)
}

# helper function that makes the function evalutaion
tf_evaluation <- function(predicted_tf, variables){
  predicted_tf_num <- predicted_tf
  for(i in variables){
    predicted_tf_num <- gsub(i, "1.0", predicted_tf_num)
  }
  tf_eval <- try({
    eval(parse(text = paste('f_test <- function() {' ,  predicted_tf_num , '}', sep='')))
    f_test()
  }, silent = TRUE)
  return(tf_eval)
}

update_tfs_in_mpr_nml <- function(run_dir, setup=NULL, parameter_list, scaling, 
                                  save_para_as_ncdf = FALSE, Training=TRUE){
  
  if(Training) {
    trainval_split <- "/Training_"
  } else {
    trainval_split <- "/Validation_"
  }
  ###
  # UPDATE NUMERIC COEFFICIENTS IN MPR.NML
  ###
  
  # define cos/sin transformation of aspect
  for (par in 1:length(parameter_list)){
    parameter_list[[par]] <- gsub("aspect_sin", 
                                  "sin( aspect*( 2.0*(3.14159265359)/360.0))",
                                  parameter_list[[par]])
    parameter_list[[par]] <- gsub("aspect_cos", 
                                  "cos(aspect*(2.0*(3.14159265359)/360.0))",
                                  parameter_list[[par]])
    parameter_list[[par]] <- gsub("lai", "lai_class", parameter_list[[par]])
    parameter_list[[par]] <- gsub("lai_class_class", "lai_class", parameter_list[[par]])
    
  }
  
  
  # predefined variables
  base <- c("dem", "aspect", "slope", "bd", "sand", "clay")
  base_plus <- c("lai", "map", "mat", "mat_range")
  mhm <- c("ThetaS", "KSat", "vGenu_n")
  lai <- c("lai_class","lai_agg")
  variables <- c(base, base_plus, mhm, lai)
  
  # get numeric coefficients of TFs
  numerics <- c()
  for(par in 1:length(parameter_list)){
    prep_tf <- prepare_tf_for_mhm(parameter_list[[par]], scaling, 
                                  variables = variables)
    numerics <- c(numerics, prep_tf$numerics)
    parameter_list[[par]] <- prep_tf$`function`
  }
  
  # make numeric vectors with parameter names
  numeric_paras <- numerics %>% unique()
  if(length(numeric_paras) > 0){
    numeric_paras_names <- paste0("GenAI_para", seq_along(numeric_paras))
    numeric_paras <- numeric_paras[order(nchar(numeric_paras), decreasing = TRUE)]
    # fill in numerics in functions
    for(par in 1:length(parameter_list)){
      for(k in seq_along(numeric_paras)){
        par_fun <- parameter_list[[par]]
        par_fun <- gsub(numeric_paras[k], numeric_paras_names[k], par_fun)
        parameter_list[[par]] <- par_fun
      }
    }
    
    # update master mpr.nml file to include all numeric coefficients
    mpr_nml <- readLines("/gpfs/data/fs71468/GenAI_para_runs/data/05_mhm/master_mhm_config/mpr.nml")
    
    # get parameter names, parameter_values
    parameter_name_line <- grep('parameter_names', mpr_nml)
    parameter_value_line <- grep('parameter_values', mpr_nml)
    parameters_end <- grep('&data_arrays', mpr_nml, fixed = TRUE) - 3
    parameter_names_initial <- mpr_nml[parameter_name_line:(parameter_value_line-1)]
    parameter_values_initial <- mpr_nml[parameter_value_line:parameters_end]
    
    # remove previous GenAI_para lines
    old_GenAI_para_nums <- parameter_names_initial[grep("GenAI_para", parameter_names_initial)]
    num_old_nums <- length(unlist(strsplit(old_GenAI_para_nums, ", ")))
    parameter_names_initial <- parameter_names_initial[!grepl("GenAI_para", parameter_names_initial)]
    all_valus <- unlist(strsplit(parameter_values_initial, ","))
    last_val <- all_valus[length(all_valus)-num_old_nums]
    last_row <- parameter_values_initial[grep(last_val, parameter_values_initial)]
    last_row_id <- which(parameter_values_initial == last_row)
    parameter_values_initial <- parameter_values_initial[1:last_row_id]
    
    
    # add new parameter names and values
    nlines <- length(parameter_names_initial)
    parameter_names_initial[nlines] <- paste0(parameter_names_initial[nlines], ",")
    parameter_names <- c(
      parameter_names_initial,
      paste0("                            ",
             paste0(paste0("'", numeric_paras_names, "'"), collapse = ", ")))
    nlines <- length(parameter_values_initial)
    parameter_values_initial[nlines] <- paste0(parameter_values_initial[nlines], ",")
    parameter_values <- c(
      parameter_values_initial,
      paste0("                            ", paste0(numeric_paras, collapse = ", ")))
    
    
    total_nr_of_numeric_paras <- length(unlist(strsplit(parameter_values, ",")))
    parameter_names[1] <- paste0("    parameter_names(1:",total_nr_of_numeric_paras, ") = ",
                                 strsplit(parameter_names[1], " = ", fixed = TRUE)[[1]][2])
    parameter_values[1] <- paste0("    parameter_values(1:",total_nr_of_numeric_paras, ") = ",
                                  strsplit(parameter_values[1], " = ", fixed = TRUE)[[1]][2])
    
    # add new parmeter names and values
    new_mpr_nml <- c(mpr_nml[1:(parameter_name_line-1)],
                     parameter_names,
                     parameter_values,
                     mpr_nml[(parameters_end+1):length(mpr_nml)])
    
    # save mpr.nml for every basin in current setup
    basins <- list.files(paste0(run_dir, trainval_split, setup))
    # adapt paths in mhm.nml for each basin and change time period
    for(basin in basins){
      # create basin specific mhm config and output folder
      basin_dir <- paste0(run_dir, trainval_split, setup, "/", basin)
      mpr_nml <- new_mpr_nml
      mpr_nml <- sapply(mpr_nml, function(x) gsub("<DATA_DIR>", data_dir, x))
      mpr_nml <- sapply(mpr_nml, function(x) gsub("<SPLIT>", "Training", x))
      mpr_nml <- sapply(mpr_nml, function(x) gsub("<BASIN>", basin, x))
      writeLines(mpr_nml, paste0(basin_dir,"/config/mpr.nml"))
    }
  }
  
  ###
  # Update TF strings in mpr.nml
  ###
  KSat <- parameter_list[["Ksat"]]
  FieldCap <- parameter_list[["FieldCap"]]
  froots_1 <- parameter_list[["fRoots_1"]]
  froots_2 <- parameter_list[["fRoots_2"]]
  ThetaS_1 <- parameter_list[["ThetaS_1"]]
  ThetaS_2 <- parameter_list[["ThetaS_2"]]
  Canopy_Intercept <- parameter_list[["L1_Max_Canopy_Intercept"]]
  
  # change lai to lai_agg in all parameter except canopy intercept
  KSat <- gsub("lai_class", "lai_agg", KSat)
  FieldCap <- gsub("lai_class", "lai_agg", FieldCap)
  froots_1 <- gsub("lai_class", "lai_agg", froots_1)
  froots_2 <- gsub("lai_class", "lai_agg", froots_2)
  ThetaS_1 <- gsub("lai_class", "lai_agg", ThetaS_1)
  ThetaS_2 <- gsub("lai_class", "lai_agg", ThetaS_2)
  
  
  # adapt till/notill versions for FieldCap and ThetaS
  till_variables <- c("sand", "clay", "bd", "KSat", "FieldCap", "ThetaS", "vGenu_n",
                      "aspect", "slope", "dem", "map", "mat", "ma*t_range", "lai_agg")
  FieldCap <- gsub("FieldCap_c1", "FC_c1", FieldCap)
  FieldCap <- gsub("FieldCap_c2", "FC_c2", FieldCap)
  FieldCap <- gsub("mat_range", "ma*t_range", FieldCap, fixed = TRUE)
  FieldCap_till <- FieldCap_notill <- FieldCap
  
  ThetaS_1 <- gsub("PTF_lower66_5_clay", "TS_c1", ThetaS_1)
  ThetaS_1 <- gsub("PTF_higher66_5_clay", "TS_c2", ThetaS_1)
  ThetaS_2 <- gsub("PTF_lower66_5_clay", "TS_c1", ThetaS_2)
  ThetaS_2 <- gsub("PTF_higher66_5_clay", "TS_c2", ThetaS_2)
  ThetaS_1 <- gsub("mat_range", "ma*t_range", ThetaS_1, fixed = TRUE)
  ThetaS_2 <- gsub("mat_range", "ma*t_range", ThetaS_2, fixed = TRUE)
  
  ThetaS_1_till <- ThetaS_1_notill <- ThetaS_1
  ThetaS_2_till <- ThetaS_2_notill <- ThetaS_2
  
  for(tvar in till_variables){
    # fieldCap
    FieldCap_till <- gsub(tvar, paste0(tvar, "_till"), FieldCap_till, fixed = TRUE)
    FieldCap_notill <- gsub(tvar, paste0(tvar, "_notill"), FieldCap_notill, fixed = TRUE)
    
    #ThetaS till (notill uses full arrays)
    if(tvar == "bd"){
      ThetaS_1_till <- gsub(tvar, paste0(tvar, "_eff_till"), ThetaS_1_till, fixed = TRUE)
      ThetaS_2_till <- gsub(tvar, paste0(tvar, "_eff_till"), ThetaS_2_till, fixed = TRUE)
    } else{
      ThetaS_1_till <- gsub(tvar, paste0(tvar, "_till"), ThetaS_1_till, fixed = TRUE)
      ThetaS_2_till <- gsub(tvar, paste0(tvar, "_till"), ThetaS_2_till, fixed = TRUE)
    }
    
    
    if(tvar == "FieldCap"){
      FieldCap_till <- gsub("FC_c1", "FieldCap_c1", FieldCap_till, fixed = TRUE)
      FieldCap_till <- gsub("FC_c2", "FieldCap_c2", FieldCap_till, fixed = TRUE)
      FieldCap_notill <- gsub("FC_c1", "FieldCap_c1", FieldCap_notill, fixed = TRUE)
      FieldCap_notill <- gsub("FC_c2", "FieldCap_c2", FieldCap_notill, fixed = TRUE)
    }
    if(tvar == "clay"){
      ThetaS_1_till <- gsub("TS_c1", "PTF_lower66_5_clay", ThetaS_1_till, fixed = TRUE)
      ThetaS_1_till <- gsub("TS_c2", "PTF_higher66_5_clay", ThetaS_1_till, fixed = TRUE)
      ThetaS_2_till <- gsub("TS_c1", "PTF_lower66_5_clay", ThetaS_2_till, fixed = TRUE)
      ThetaS_2_till <- gsub("TS_c2", "PTF_higher66_5_clay", ThetaS_2_till, fixed = TRUE)
      
      ThetaS_1_notill <- gsub("TS_c1", "PTF_lower66_5_clay", ThetaS_1_notill, fixed = TRUE)
      ThetaS_1_notill <- gsub("TS_c2", "PTF_higher66_5_clay", ThetaS_1_notill, fixed = TRUE)
      ThetaS_2_notill <- gsub("TS_c1", "PTF_lower66_5_clay", ThetaS_2_notill, fixed = TRUE)
      ThetaS_2_notill <- gsub("TS_c2", "PTF_higher66_5_clay", ThetaS_2_notill, fixed = TRUE)
    }
  }
  
  # fix mat
  FieldCap_till  <- gsub("ma*t_range", "mat_range", FieldCap_till, fixed = TRUE)
  FieldCap_notill  <- gsub("ma*t_range", "mat_range", FieldCap_notill, fixed = TRUE)
  ThetaS_1_till <- gsub("ma*t_range", "mat_range", ThetaS_1_till, fixed = TRUE)
  ThetaS_1_notill <- gsub("ma*t_range", "mat_range", ThetaS_1_notill, fixed = TRUE)
  ThetaS_2_till <- gsub("ma*t_range", "mat_range", ThetaS_2_till, fixed = TRUE)
  ThetaS_2_notill <- gsub("ma*t_range", "mat_range", ThetaS_2_notill, fixed = TRUE)
  
  # change KSat variables to _horizon if necessary
  for(var in c("aspect", "slope", "dem", "mat", "mat_range", "map", "lai_agg")){
    KSat <- gsub(var, paste0(var, "_horizon"), KSat)
    ThetaS_1_notill <- gsub(var, paste0(var, "_horizon"), ThetaS_1_notill)
    ThetaS_2_notill <- gsub(var, paste0(var, "_horizon"), ThetaS_2_notill)
    KSat <- gsub("mat_horizon_range", "mat_range", KSat)
    ThetaS_1_notill <- gsub("mat_horizon_range", "mat_range", ThetaS_1_notill)
    ThetaS_2_notill <- gsub("mat_horizon_range", "mat_range", ThetaS_2_notill)
  }
  
  for(var in c("clay", "sand", "bd", "aspect", "slope", "dem", "mat", "mat_range", "map", "lai_agg")){
    froots_1 <- gsub(var, paste0(var, "_land_cover"), froots_1)
    froots_1 <- gsub("mat_land_cover_range", "mat_range", froots_1)
    
    froots_2 <- gsub(var, paste0(var, "_land_cover"), froots_2)
    froots_2 <- gsub("mat_land_cover_range", "mat_range", froots_2)
    
  }
  
  
  # creat full TFs from parts
  fRoots <- paste0("where (land_cover_horizon > 0.5 .and. land_cover_horizon < 1.5)(",
                   froots_1,
                   ") else where (land_cover_horizon < 2.5)(1.0 - rootFractionCoefficient_impervious ** (z_upper_bound * 100.0)) - (1.0 - rootFractionCoefficient_impervious ** (z_lower_bound * 100.0)) else where (land_cover_horizon < 3.5) (",
                   froots_2, ")")
  ThetaS_till <- paste0("where (sand_till < vGenu_tresh) (",
                        ThetaS_1_till, 
                        ") else where (sand_till >= vGenu_tresh) (",
                        ThetaS_2_till, ")")
  ThetaS_notill <- paste0("where (sand < vGenu_tresh) (",
                          ThetaS_1_notill, 
                          ") else where (sand >= vGenu_tresh) (",
                          ThetaS_2_notill, ")")
  # detect inputs
  all_inputs <- c("KSat", "KSat_till", "KSat_notill", 
                  "ThetaS_till", "ThetaS_notill", 
                  "vGenu_n_till", "vGenu_n_notill",
                  "aspect", "slope", "dem",
                  "aspect_till", "slope_till", "dem_till",
                  "aspect_notill", "slope_notill", "dem_notill",
                  "aspect_horizon", "slope_horizon", "dem_horizon", 
                  "sand", "clay", "bd",
                  "sand_till", "sand_notill", "clay_till", "clay_notill", "bd_till", "bd_notill",
                  "map", "mat", "mat_range",
                  "map_till", "mat_till", "mat_range_till",
                  "map_notill", "mat_notill", "mat_range_notill",
                  "map_horizon", "mat_horizon", "mat_range_horizon", 
                  "lai_agg", "lai_agg_horizon", "lai_agg_till", "lai_agg_notill",
                  "lai_class", "bd_eff_till", "land_cover_horizon", "z_lower_bound", "z_upper_bound",
                  "clay_land_cover", "sand_land_cover", "bd_land_cover",
                  "map_land_cover", "mat_land_cover", "mat_range_land_cover", "lai_agg_land_cover",
                  "aspect_land_cover", "slope_land_cover", "dem_land_cover")
  
  tmp_vars <- unlist(strsplit(gsub(" ", "", KSat), "(?=[+-/*><=)()])", perl = TRUE))
  KSat_inputs <- paste0(unique(tmp_vars[tmp_vars %in% all_inputs]), collapse=", ")
  
  tmp_vars <- unlist(strsplit(gsub(" ", "", FieldCap_till), "(?=[+-/*><=)()])", perl = TRUE))
  FieldCap_till_inputs <- paste0(unique(tmp_vars[tmp_vars %in% all_inputs]), collapse=", ")
  tmp_vars <- unlist(strsplit(gsub(" ", "", FieldCap_notill), "(?=[+-/*><=)()])", perl = TRUE))
  FieldCap_notill_inputs <- paste0(unique(tmp_vars[tmp_vars %in% all_inputs]), collapse=", ")
  
  tmp_vars <- unlist(strsplit(gsub(" ", "", fRoots), "(?=[+-/*><=)()])", perl = TRUE))
  froots_inputs <- paste0(unique(tmp_vars[tmp_vars %in% all_inputs]), collapse=", ")
  
  tmp_vars <- unlist(strsplit(gsub(" ", "", ThetaS_till), "(?=[+-/*><=)()])", perl = TRUE))
  ThetaS_till_inputs <- paste0(unique(tmp_vars[tmp_vars %in% all_inputs]), collapse=", ")
  tmp_vars <- unlist(strsplit(gsub(" ", "", ThetaS_notill), "(?=[+-/*><=)()])", perl = TRUE))
  ThetaS_notill_inputs <- paste0(unique(tmp_vars[tmp_vars %in% all_inputs]), collapse=", ")
  
  tmp_vars <- unlist(strsplit(gsub(" ", "", Canopy_Intercept), "(?=[+-/*><=)()])", perl = TRUE))
  Canopy_Intercept_inputs <- paste0(unique(tmp_vars[tmp_vars %in% all_inputs]), collapse=", ")
  
  run_tfs <- data.frame(parameters = c("KSat", "FieldCap_till", "FieldCap_notill", "ThetaS_till", 
                                       "ThetaS_notill", "fRoots_temp", "L1_Max_Canopy_Intercept"),
                        TFs = c(KSat, FieldCap_till, FieldCap_notill, 
                                ThetaS_till, ThetaS_notill, fRoots, Canopy_Intercept),
                        inputs = c(KSat_inputs, FieldCap_till_inputs, FieldCap_notill_inputs,
                                   ThetaS_till_inputs, ThetaS_notill_inputs, froots_inputs,
                                   Canopy_Intercept_inputs))
  run_tfs$TFs <- gsub("(", "( ", run_tfs$TFs, fixed = TRUE)
  run_tfs$TFs <- gsub("(  ", "( ", run_tfs$TFs, fixed = TRUE)
  if(is.null(setup)){
    write.csv(run_tfs, paste0(run_dir, "/current_tfs.csv"), row.names = FALSE)
    
    # write relevant mpr parameters to ncdf?
    save_para <- ifelse(save_para_as_ncdf, "True", "False")
    
    # update all mpr files with python script
    system(paste0("/home/fs71468/mfeigl/miniconda3/envs/GenAI_para_mhm_py3/bin/python3",
                  " /gpfs/data/fs71468/GenAI_para_runs/code/02_utility_functions/08_update_tfs_in_mpr_nml.py -f ",
                  run_dir, trainval_split, " -t ",
                  run_dir, "/current_tfs.csv -s ", save_para))
  } else {
    write.csv(run_tfs, paste0(run_dir, "/current_tfs", setup, ".csv"), row.names = FALSE)
    
    # write relevant mpr parameters to ncdf?
    save_para <- ifelse(save_para_as_ncdf, "True", "False")
    
    # update all mpr files with python script
    
    sys_out <- robust.system(
      paste0("/home/fs71468/mfeigl/miniconda3/envs/GenAI_para_mhm_py3/bin/python3",
             " /gpfs/data/fs71468/GenAI_para_runs/code/02_utility_functions/08_update_tfs_in_mpr_nml.py -f ",
             run_dir, trainval_split, setup, " -t ",
             run_dir, "/current_tfs", setup, ".csv -s ", save_para))
    if(length(grep("Error", sys_out$stderr)) > 0){
      stop("Failed TF update in mHM source code")
    }
  }
}


update_numeric_parameters <- function(numeric_parameters, 
                                      run_dir,
                                      setup=NULL,
                                      parameter_nml=parameter_nml, 
                                      ind_of_num_paras=ind_of_num_paras,
                                      Training=TRUE){
  
  if(Training) {
    trainval_split <- "/Training"
  } else {
    trainval_split <- "/Validation"
  }
  
  for(i_numerics in seq_along(numeric_parameters)){
    i_para <- ind_of_num_paras[i_numerics]
    parameter_nml$mhm_parameters[[i_para]][3] <- numeric_parameters[i_numerics]
  }
  if(is.null(setup)){
    all_files <- paste0(list.files(paste0(run_dir, trainval_split), full.names = TRUE),
                        "/config/mhm_parameter.nml")
  } else {
    all_files <- paste0(list.files(paste0(run_dir, trainval_split, "_", setup), full.names = TRUE),
                        "/config/mhm_parameter.nml")
  }
  tmp_write_fun <- function(file){
    write_mhm_nml(parameter_nml, file)
  }
  no_return <- parallel::mclapply(all_files, tmp_write_fun, mc.cores = length(all_files)/3)
}


SPAEF <- function(observations, predictions){
  spaef_try <- try({
    
    # correlation
    alpha <- cor(predictions, observations)
    
    # coefficient of variation
    cv_obs <- sd(observations)/mean(observations)
    cv_pred <- sd(predictions)/mean(predictions);
    beta <- cv_pred/cv_obs;
    
    # histogram distance
    observations <- (observations-mean(observations))/sd(observations)
    predictions <- (predictions-mean(predictions))/sd(predictions)
    bins <- floor(sqrt(length(observations)))
    suppressWarnings({
      h1 <- hist(observations, breaks = bins, freq = TRUE, plot = FALSE)
      h2 <- hist(predictions, breaks = bins, freq = TRUE, plot = FALSE) 
    })
    a <- pracma::histc(observations, h1$breaks)
    b <- pracma::histc(predictions, h1$breaks)
    c <- cbind(a$cnt, b$cnt)
    d <- pmin(c[,1],c[,2])
    overlap <- sum(d)
    gamma <- overlap/sum(a$cnt)
    spaef <- 1 - sqrt((alpha-1)^2 + (beta-1)^2 + (gamma-1)^2)
  }, silent = TRUE)
  if(class(spaef_try) == "try-error") spaef <- NA
  return(spaef)
}

basins_quality_criteria <- function(run_dir, setup, data_dir, Training=TRUE, ET=TRUE){
  # computes NSE, log NSE, KGE and SPAEF for provided basins
  # input: chr, basin directories
  # output: data.frame, quality criteria for each basins as data frame
  if(is.null(setup)){
    setup <- ""
  } else {
    setup <- paste0("_", setup)
  }
  if(Training) {
    trainval_split <- "/Training"
  } else {
    trainval_split <- "/Validation"
  }
  # get ET and Q results and compute quality criteria
  basins <- list.files(paste0(run_dir, trainval_split, setup))
  basin_qc <- data.frame("Basin" = basins, 
                         "KGE" = NA, 
                         "NSE" = NA, 
                         "lNSE" = NA, 
                         "SPAEF" = NA)
  qc_try <- try({
    
    for(basin in basins){
      
      if(ET){
        # ET mhm
        mhm_et_data <- nc_open(paste0(run_dir, trainval_split, setup, "/", basin,
                                      "/output/mHM_Fluxes_States.nc"))
        mhm_et_nc <-  ncvar_get(mhm_et_data, "aET")
        # get only mean of summer months May-Sep and average over years
        years <- dim(mhm_et_nc)[3] / 12
        relevant_months <- 5:9
        mhm_monthly_means <- array(NA, dim = dim(mhm_et_nc[, , relevant_months]))
        for(month in seq_along(relevant_months)){
          months <- relevant_months[month]
          for(i in 2:years) months <- c(months, relevant_months[month] +(12*(i-1)))
          mhm_monthly_means[, , month] <- apply(mhm_et_nc[, , months], 1:2, mean)
        }
        
        # ET Alexi - compute means for each month (i.e. average over years)
        et_data <- nc_open(paste0(data_dir, "06_training_validation_data/", trainval_split,
                                  "/", basin, "/ET/et.nc"))
        et_nc <- ncvar_get(et_data, "__xarray_dataarray_variable__")
        
        relevant_months <- 1:5
        alexi_monthly_means <- array(NA, dim = dim(et_nc[, , relevant_months]))
        for(month in relevant_months){
          months <- relevant_months[month]
          for(i in 2:years) months <- c(months, month +(5*(i-1)))
          alexi_monthly_means[, , month] <- apply(et_nc[, , months], 1:2, mean)
        }
        
        
        # SPAEF
        spaef_months <- numeric(dim(mhm_monthly_means)[3])
        for(month in 1:length(spaef_months)){
          # get ET values of specific month as vector
          sim_et <- as.numeric(mhm_monthly_means[, , month])
          obs_et <- as.numeric(alexi_monthly_means[, , month])
          # remove values that are NA in the simulated er
          obs_et <- obs_et[!is.na(sim_et)]
          sim_et <- sim_et[!is.na(sim_et)]
          # remove values that are NAN in the observed et
          sim_et <- sim_et[!is.nan(obs_et)]
          obs_et <- obs_et[!is.nan(obs_et)]
          # compute SPAEFF
          spaef_months[month] <- SPAEF(obs_et, sim_et)
        }
        nc_close(mhm_et_data)
        nc_close(et_data)
        spaef <- mean(spaef_months)
      } else {
        spaef <- NA
      }
      # Discharge
      basin_q <- read.table(paste0(run_dir, trainval_split, setup, "/", basin, 
                                   "/output/daily_discharge.out"), header = TRUE)
      basin_q <- basin_q[, c(5, 6)]
      names(basin_q) <- c("Qobs", "Qsim")
      basin_q[which(basin_q < 0, arr.ind = TRUE)] <- NA
      # Discharge criteria
      kge <- KGE(basin_q$Qsim, basin_q$Qobs)
      nse <- NSE(basin_q$Qsim, basin_q$Qobs)
      log_basin_q <- suppressWarnings(log(basin_q))
      inf_rows <- c(which(is.infinite(log_basin_q$Qsim)), which(is.infinite(log_basin_q$Qobs)))
      log_basin_q[inf_rows, ] <- NA 
      lnse <- suppressWarnings(NSE(log_basin_q$Qsim, log_basin_q$Qobs))
      # write results
      basin_qc[basin_qc$Basin == tail(strsplit(basin, "/")[[1]], 1), 
               c("KGE", "NSE", "lNSE", "SPAEF")] <- c(kge, nse, lnse, spaef)
    }
  })
  return(basin_qc)
}

# computes loss given the quality criteria and the mode
compute_loss <- function(basin_qc, mode, point_tf = NULL){
  
  loss_try <- try({
    # NSE, log NSE
    if(mode == "Q"){
      multi_obj_loss <- sign(basin_qc$NSE) * 1/2*abs(basin_qc$NSE)^6 +
        sign(basin_qc$lNSE) * 1/2*abs(basin_qc$lNSE)^6
      size_loss <- 0
    }
    
    # NSE, log NSE, SPAEF
    if(mode == "QET"){
      multi_obj_loss <- sign(basin_qc$NSE) * 1/3*abs(basin_qc$NSE)^6 +
        sign(basin_qc$lNSE) * 1/3*abs(basin_qc$lNSE)^6 +
        sign(basin_qc$SPAEF) * 1/3*abs(basin_qc$SPAEF)^6
      size_loss <- 0
    }
    
    # NSE, log NSE, TF size
    if(mode == "Q_TF"){
      multi_obj_loss <- sign(basin_qc$NSE) * 1/3*abs(basin_qc$NSE)^6 +
        sign(basin_qc$lNSE) * 1/3*abs(basin_qc$lNSE)^6
      functions_splitted <- lapply(point_tf, .function_splitter)
      size_loss <- length(unlist(functions_splitted)) * 0.001/2
    }
    
    # NSE, log NSE, SPAEF, TF size
    if(mode == "QET_TF"){
      multi_obj_loss <- sign(basin_qc$NSE) * 1/3*abs(basin_qc$NSE)^6 +
        sign(basin_qc$lNSE) * 1/3*abs(basin_qc$lNSE)^6 +
        sign(basin_qc$SPAEF) * 1/3*abs(basin_qc$SPAEF)^6
      functions_splitted <- lapply(point_tf, .function_splitter)
      size_loss <- length(unlist(functions_splitted)) * 0.001/2
    }
    
    # Loss
    multi_obj_loss <- ifelse(multi_obj_loss < 0,
                             -1*(multi_obj_loss*-1)^(1/6),
                             multi_obj_loss^(1/6))
    wmulti_obj_loss <- mean(x = multi_obj_loss)
    loss <- wmulti_obj_loss - size_loss
    
  })
  
  # catch errors in loss computation
  if(class(loss_try) == "try-error"){
    loss <- NA
    domain_wloss <- matrix(NA, ncol = nrow(basin_qc))
    colnames(domain_wloss) <- basin_qc$Basin
  }
  
  # Loss per basin
  domain_wloss <- matrix(multi_obj_loss, ncol = nrow(basin_qc))
  colnames(domain_wloss) <- basin_qc$Basin
  
  return(list(loss, domain_wloss))
}

# Runs mHM on multiple nodes
run_mhm <- function(run_dir, setup = 1, mhm_preparation = TRUE, Training=TRUE){
  
  if(Training) {
    trainval_split <- "/Training_"
    run_script_folder <- "training"
  } else {
    trainval_split <- "/Validation_"
    run_script_folder <- "Validation"
  }
  
  # Run python script and compile mhm
  if(mhm_preparation){
    system(paste0("bash 05_run_mhm_preparation_", setup, ".sh & wait"))
  }
  
  # remove previous flux_states if they exist
  basins <- list.files(paste0(run_dir, trainval_split, setup))
  for (basin in basins){
    unlink(paste0(run_dir, trainval_split, setup, "/", basin, "/output/mHM_Fluxes_States.nc"))
  }
  system(paste0("export LD_LIBRARY_PATH=$LIBRARY_PATH
  cd ", run_script_folder, "_run_scripts_", setup, "
                  for j in *
                  do
                  bash $j &
                  done
                  wait"))
}


track_results <- function(setup, iter, result_dir, point, basin_qc, loss_list, 
                          num_para_bounds, ind_of_num_paras,
                          point_tf = NULL, TF_names, only_numerics, latent_dim){
  # track results in 5 files
  # TFs_tracker:    iteration, TF1, TF2, ...
  # losses_tracker: iteration, loss, mean NSE, mean log NSE, mean KGE, mean SPAEF
  # basin_losses_tracker: iteration, basin_qc data frame
  # parameters_tracker: iteration, point with correct parameter names
  # best_TFs_tracker: iteration, best TFs & numeric parameters
  
  point <- as.numeric(point)
  # initialization
  if(!file.exists(paste0(result_dir, "/parameters_tracker_",setup, ".csv"))){
    TFs_tracker <- NULL
    losses_tracker <- NULL
    basin_qc_tracker <- NULL
    parameters_tracker <- NULL
    best_TFs_tracker <- NULL
  } else {
    if(!only_numerics) TFs_tracker <- read.csv(paste0(result_dir, "/TFs_tracker_",setup, ".csv"))
    losses_tracker <- read.csv(paste0(result_dir, "/losses_tracker_",setup, ".csv"))
    basin_qc_tracker <- feather::read_feather(paste0(result_dir, "/basin_qc_tracker_",setup, ".feather"))
    parameters_tracker <- read.csv(paste0(result_dir, "/parameters_tracker_",setup, ".csv"))
    best_TFs_tracker <- read.csv(paste0(result_dir, "/best_TFs_tracker_",setup, ".csv"))
  }
  
  # timestamp
  timestamp <- as.character(Sys.time())
  
  # TFs
  if(!only_numerics){
    TFs_tracker_add <- data.frame(iteration = iter,
                                  complex = setup,
                                  time = timestamp,
                                  matrix(point_tf, ncol = length(point_tf)))
    colnames(TFs_tracker_add)[-c(1:3)] <- TF_names
  }
  # losses
  losses_tracker_add <- data.frame(iteration = iter,
                                   complex = setup,
                                   time = timestamp,
                                   loss = loss_list[[1]],
                                   mean_NSE = mean(basin_qc$NSE),
                                   mean_lNSE = mean(basin_qc$lNSE),
                                   mean_KGE = mean(basin_qc$KGE),
                                   mean_SPAEF = mean(basin_qc$SPAEF))
  
  
  # basine wise quality criteria
  basin_qc_tracker_add <- cbind(iteration = iter,
                                complex = setup,
                                time = timestamp,
                                basin_qc)
  
  # parameters -> numeric vector that is optimized
  parameters_tracker_add <- data.frame(iteration = iter,
                                       complex = setup,
                                       time = timestamp,
                                       matrix(point, ncol = length(point)))
  if(only_numerics){
    names(parameters_tracker_add)[-c(1:3)] <- row.names(num_para_bounds)
  } else {
    names(parameters_tracker_add)[-c(1:3)] <- c(paste0("GenAI_para", 1:(latent_dim*length(point_tf))),
                                                row.names(num_para_bounds))
  }
  # best TF results
  best_loss <- ifelse(is.null(best_TFs_tracker), -9999, max(best_TFs_tracker$loss, na.rm=TRUE))
  best_loss <- ifelse(best_loss == -Inf, -9999, best_loss)
  loss_proxy <- ifelse(is.na(loss_list[[1]]), -9998, loss_list[[1]])
  if(loss_proxy > best_loss){
    if(only_numerics){
      best_TFs_tracker_add <- data.frame(iteration = iter,
                                         complex = setup, 
                                         time = timestamp,
                                         loss = loss_list[[1]])
      best_TFs_tracker_add <- cbind(best_TFs_tracker_add, 
                                    parameters_tracker_add[, -c(1:3)])
    } else {
      best_TFs_tracker_add <- data.frame(iteration = iter,
                                         complex = setup, 
                                         time = timestamp,
                                         loss = loss_list[[1]],
                                         matrix(point_tf, ncol = length(point_tf)))
      names(best_TFs_tracker_add)[-c(1:4)] <- TF_names
      best_TFs_tracker_add <- cbind(best_TFs_tracker_add, 
                                    parameters_tracker_add[, -c(1:3)])
    }
  }
  
  # set names in case the csv name formats changed
  if(!is.null(best_TFs_tracker)){
    if(!only_numerics) names(TFs_tracker_add) <- names(TFs_tracker)
    names(losses_tracker_add) <- names(losses_tracker)
    names(parameters_tracker_add) <- names(parameters_tracker)
    names(basin_qc_tracker_add) <- names(basin_qc_tracker)
    if(exists("best_TFs_tracker_add")) names(best_TFs_tracker_add) <- names(best_TFs_tracker)
  }
  
  if(!only_numerics){
    write.csv(rbind(TFs_tracker, TFs_tracker_add), 
              paste0(result_dir, "/TFs_tracker_",setup, ".csv"), row.names = FALSE)
  }
  write.csv(rbind(losses_tracker, losses_tracker_add), 
            paste0(result_dir, "/losses_tracker_",setup, ".csv"), row.names = FALSE)
  feather::write_feather(rbind(basin_qc_tracker, basin_qc_tracker_add), 
                         paste0(result_dir, "/basin_qc_tracker_",setup, ".feather"))
  parameters_tracker <- rbind(parameters_tracker, parameters_tracker_add)
  write.csv(parameters_tracker, 
            paste0(result_dir, "/parameters_tracker_", setup, ".csv"), row.names = FALSE)
  if(exists("best_TFs_tracker_add")){
    best_TFs_tracker <- rbind(best_TFs_tracker, best_TFs_tracker_add)
  }
  write.csv(best_TFs_tracker, 
            paste0(result_dir, "/best_TFs_tracker_",setup, ".csv"), row.names = FALSE)
}


objective_function_definer <- function(run_dir, 
                                       data_dir,
                                       only_numerics,
                                       parameter_nml,
                                       ind_of_num_paras,
                                       loss_mode, 
                                       latent_dim = NULL){
  
  numeric_optim_obj <- function(point, setup, iter){
    
    # Make a counter for the optimization
    cat(iter, "iteration -", setup, "complex\n")
    
    # update numeric parameters
    update_numeric_parameters(point, run_dir, setup, parameter_nml, ind_of_num_paras)
    
    # Run mHM
    start_eval <- Sys.time()
    run_mhm(run_dir, setup, mhm_preparation = FALSE)
    end_eval <- Sys.time()
    cat("Evaluation took:\n")
    print(end_eval - start_eval)
    
    # compute quality criteria and loss
    basin_qc <- basins_quality_criteria(run_dir, setup, data_dir)
    loss_list <- compute_loss(basin_qc, mode = loss_mode)
    
    # track results
    track_results(setup, iter, result_dir, point, basin_qc, loss_list,
                  num_para_bounds, ind_of_num_paras,
                  point_tf = NULL, TF_names, only_numerics = TRUE, latent_dim=latent_dim)
    
    # return loss
    return(loss_list[[1]])
  }
  
  tf_optim_obj <- function(point, setup, iter){
    
    # Make a counter for the optimization
    cat("Complex", setup, "/ Iteration",iter, "\n")
    # Update TF and numeric coefficients
    vae_pred_path <- paste0(run_dir, "/current_vae_predictions/")
    dir.create(vae_pred_path, showWarnings = FALSE)
    # store vae inputs for python script
    start_time <- Sys.time()
    parameter_list <- list()
    if(class(point) != "data.frame"){
      point <- as.data.frame(matrix(as.numeric(point), ncol=length(point)))
    }
    
    # generate TF using python scripts running in background
    loop_over_tfs <- names(vae_list)
    for(TF in loop_over_tfs){
      parameter_list[TF] <- NA
      tf_generate_count <- 0
      retry_tf_generation <- TRUE
      while(retry_tf_generation){
        tf_generate_count <- tf_generate_count + 1
        # define and write vae inputs
        range <- switch(TF,
                        "Ksat" = 1:30, "FieldCap" = 31:60,"ThetaS_1" = 61:90,
                        "ThetaS_2" = 91:120,"fRoots_1" = 121:150,
                        "fRoots_2" = 151:180,"L1_Max_Canopy_Intercept" = 181:210)
        vae <- vae_list[[TF]]
        status_file <- paste0(vae_pred_path, "vae", vae, "_setup", setup, "_start.txt")
        write.csv(point[range], 
                  paste0(vae_pred_path, "vae", vae, "_setup", setup, "_VAEinput.csv"), 
                  row.names = FALSE)
        writeLines("start", status_file) # start key for vae script in background
        
        wait_for_vae <- TRUE
        start_while <- Sys.time()
        while(wait_for_vae){
          try({if(readLines(status_file) == "end"){
            wait_for_vae <- FALSE}}, silent = TRUE)
          
          if(difftime(Sys.time(), start_while, units = "secs")[[1]] > 50){
            cat("restart vae...\n")
            system("killall python3")
            start_GenAI_para_vae(run_dir, code_dir, vae_list, ncomplexes)
            writeLines("start", status_file)
            start_while <- Sys.time()
          }
        }
        # compute TF from softmax prediction
        softmax_pred <- read.csv(
          paste0(run_dir, "/current_vae_predictions/vae", vae, "_setup", setup, "_softmax.csv"))
        parameter_list[TF] <- generate_function_from_softmax(vae_list[[TF]], 
                                                             softmax_pred, 
                                                             data_dir, variables)
        # retry TF generation if NA generated, stop if valid TF generated or after 10 tries
        if(!is.na(parameter_list[TF]) | tf_generate_count > 30) retry_tf_generation <- FALSE
      }
      # scale parameters if TFs map outside their ranges
      if(TF %in% names(scale_list)){
        quant_pred <- read.csv(
          paste0(run_dir, "/current_vae_predictions/vae", vae, "_setup", setup, "_quantiles.csv"))
        if(quant_pred[, 9] > scale_list[[TF]][2]){
          parameter_list[TF] <- paste0("( ", parameter_list[TF], " ) / ", scale_factor)
        }
      }
    }
    names(parameter_list) <- names(vae_list)
    
    # fix names
    point_tf <- do.call(c, parameter_list)
    TF_names <- names(point_tf)
    names(point_tf) <- NULL
    
    # if TFs without variables are generated --> NA
    for(l in seq_along(point_tf)){
      if(sum(sapply(variables, function(x) grepl(x, point_tf[l]))) == 0) point_tf[l] <- NA
      if(grepl("NA", point_tf[l])) point_tf[l] <- NA
    }
    
    # return NA loss in case that invalid TFs are generated
    if(sum(is.na(point_tf)) > 0){
      basins <- list.files(paste0(run_dir, "/Training_", setup))
      basin_qc <- data.frame("Basin" = basins, 
                             "KGE" = NA, 
                             "NSE" = NA, 
                             "lNSE" = NA, 
                             "SPAEF" = NA)
      loss_list <- compute_loss(basin_qc, loss_mode, point_tf)
      track_results(setup, iter, result_dir, point, basin_qc, loss_list, 
                    num_para_bounds, ind_of_num_paras,
                    point_tf, TF_names, only_numerics = FALSE, latent_dim=latent_dim)
      system(paste0("chmod -R 777 ", run_dir))
      system(paste0("chmod -R 777 ", result_dir))
      return(loss_list[[1]])
    }
    
    # if all TFs are valid run mHM and compute losses
    update_tfs_in_mpr_nml(run_dir, setup, parameter_list, scaling = FALSE, 
                          save_para_as_ncdf = FALSE)
    # numeric para update
    numeric_para_ind <- (length(parameter_list)*latent_dim + 1):length(point)
    update_numeric_parameters(point[numeric_para_ind], run_dir, setup, 
                              parameter_nml, ind_of_num_paras)
    # Run mHM
    run_mhm(run_dir, setup, mhm_preparation = TRUE)
    
    # compute quality criteria and loss
    basin_qc <- basins_quality_criteria(run_dir, setup, data_dir)
    loss_list <- compute_loss(basin_qc, loss_mode, point_tf)
    
    # track results
    track_results(setup, iter, result_dir, point, basin_qc, loss_list, 
                  num_para_bounds, ind_of_num_paras,
                  point_tf, TF_names, only_numerics = FALSE, latent_dim=latent_dim)
    
    end_time <- Sys.time()
    # print status
    for(tf in TF_names){
      cat(paste0(tf, ":"), parameter_list[[tf]], "\n")
    }
    cat("Loss:", loss_list[[1]], "\n")
    cat("Run time:") 
    print(end_time - start_time)
    
    return(loss_list[[1]])
  }
  
  
  if(only_numerics) return(numeric_optim_obj)
  if(!only_numerics) return(tf_optim_obj)
}

sample_initial_popoulation <- function(run_dir, experiment, obj, ncomplexes, NDIM, 
                                       vae_list, latent_dim, standard_parameters){
  # Select initial population in case of experiment=[1,2] and no previous iterations
  if(experiment %in% 1:2 & obj$iterations == 0 & !("initial_population.csv" %in% list.files(run_dir))){
    cat("Sample initial population\n")
    nPOINTS <- ncomplexes*2*NDIM + 1
    ini_pop_size <- nPOINTS*50 # sample additional points in case NAs are produced
    POPULATION <- data.frame(matrix(NA, nrow = ini_pop_size, ncol = NDIM))
    for(k in 1:length(vae_list)){
      parameter <- names(vae_list)[k]
      vae <- vae_list[[parameter]]
      all_pop <- read.csv(paste0(data_dir, "07_VAE_data/04_vae_data/vae_", vae, 
                                 "/parameter_specific/","vae", vae, "_", parameter, 
                                 "_samples.csv"))
      sample_ids <- sample(nrow(all_pop), size = ini_pop_size)
      POPULATION[, (1:latent_dim)+((k-1)*latent_dim)] <- all_pop[sample_ids, 
                                                                 paste0("LS", 1:latent_dim)]
    }
    # add numeric parameter initial values
    non_tf_columns <- (length(vae_list)*latent_dim + 1): NDIM
    POPULATION[1, non_tf_columns] <- standard_parameters
    for (i in seq_along(non_tf_columns)){
      POPULATION[-1, non_tf_columns[i]] <- runif(ini_pop_size-1, 
                                                 num_para_bounds[i, 1], 
                                                 num_para_bounds[i, 2])
    }
    write.csv(POPULATION, paste0(run_dir, "/initial_population.csv"), row.names = FALSE)
  }
  
  if(obj$iterations == 0 & "initial_population.csv" %in% list.files(run_dir)){
    POPULATION <- read.csv(paste0(run_dir, "/initial_population.csv"))
  } else {
    POPULATION <- NULL
  }
  return(POPULATION)
}












# Function to estimate SCE runtime
sce_iteration_runtime <- function(nCOMPLEXES, NDIM, eval_time, iterations = 1){
  nPOINTS_COMPLEX <- 2 * NDIM + 1
  nPOINTS_SIMPLEX <- NDIM+1
  nPOINTS <- nCOMPLEXES * nPOINTS_COMPLEX
  CCEITER <- 2 * NDIM + 1
  
  init_pop_eval <- (nCOMPLEXES * nPOINTS_COMPLEX)
  max_inner_loop_evals <- nCOMPLEXES * CCEITER * 3
  min_inner_loop_evals <- nCOMPLEXES * CCEITER * 1
  max_evals <- init_pop_eval + max_inner_loop_evals
  min_evals <- init_pop_eval + min_inner_loop_evals
  cat("Evaluations:", min_evals* iterations, " - ", max_evals* iterations, "\n")
  cat("Run time:", min_evals*eval_time/60 *iterations, " - ", 
      max_evals*eval_time/60 * iterations,  " hours\n")
  cat("Run time:", min_evals*eval_time/60/24 *iterations, " - ", 
      max_evals*eval_time/60/24 * iterations,  " days\n")
}
# sce_iteration_runtime(1, 39, 27/60, 20)
# sce_iteration_runtime(5, 1, 1, 10)


# SCE
# The function is adapted from the hydromad package (https://github.com/floybix/hydromad/)
# this version allows saving and loading of states of the optimizer for working on time 
# limited cluster ressources

# Duan et al. standard parameters
# n ... dimension of the problem
# m = 2n+1 , number of points in each subcomplex -> cce.iter
# q = n + 1 number of points in each subcomplex
# p >= 2, the number of complexes
# alpha = 1, the number of original parents that will be replaced
# beta = 2n+1, the number of offspring that should be generated

SCEoptim <- function(FUN, par,
                     lower = -Inf, upper = Inf,
                     control = list(), state_dir = NULL, 
                     restart_after_n_iter, ini_pop = NULL, ...) {
  cat("Starting SCE Optimization.\n")
  # Initial setup (always) 
  FUN <- match.fun(FUN)
  stopifnot(is.numeric(par))
  stopifnot(length(par) > 0)
  stopifnot(is.numeric(lower))
  stopifnot(is.numeric(upper))
  ## allow `lower` or `upper` to apply to all parameters
  if (length(lower) == 1)
    lower <- rep(lower, length = length(par))
  if (length(upper) == 1)
    upper <- rep(upper, length = length(par))
  stopifnot(length(lower) == length(par))
  stopifnot(length(upper) == length(par))
  ## determine number of variables to be optimized
  NDIM <- length(par)
  ## update default options with supplied options
  stopifnot(is.list(control))
  control <- modifyList(sceDefaults(), control)
  isValid <- names(control) %in% names(sceDefaults())
  if (any(!isValid))
    stop("unrecognised options: ",
         toString(names(control)[!isValid]))
  returnpop <- control$returnpop
  trace <- control$trace
  nCOMPLEXES <- control$ncomplex
  CCEITER <- control$cce.iter
  MAXIT <- control$maxit
  MAXEVAL <- control$maxeval
  fnscale <- control$fnscale
  GenAI_para_dims <- control$GenAI_para_dims
  time_limit <- control$time_limit
  ## recommended number of CCE steps in Duan et al 1994:
  if (is.na(CCEITER))
    CCEITER <- 2 * NDIM + 1
  if (is.finite(MAXEVAL)) {
    ## upper bound on number of iterations to reach MAXEVAL
    MAXIT <- min(MAXIT, ceiling(MAXEVAL / (nCOMPLEXES * CCEITER)))
  }
  ## define number of points in each complex
  nPOINTS_COMPLEX <- control$nPOINTS_COMPLEX
  ## define number of points in each simplex
  nPOINTS_SIMPLEX <- 40 #NDIM+1
  ## define total number of points
  nPOINTS <- nCOMPLEXES * nPOINTS_COMPLEX
  
  # function definitions
  costFunction <- function(FUN, setup, iter, par, fnscale) {
    ## check lower and upper bounds
    i <- which(par < lower)
    if (any(i)) {
      cat("\nParameters", paste0("[", paste0(i, collapse = ","), "]"), "below bound\n")
      i <- i[1]
      return( 1e12 + (lower[i] - par[i]) * 1e6 )
    }
    i <- which(par > upper)
    if (any(i)) {
      cat("\nParameters", paste0("[", paste0(i, collapse = ","), "]"), "below bound\n")
      i <- i[1]
      return( 1e12 + (par[i] - upper[i]) * 1e6   )
    }
    funevals <<- funevals + 1
    fun_try <- try({
      result <- FUN(par, setup, iter) * fnscale
    })
    if(class(fun_try) == "try-error"){result <- 1e12}
    if(is.na(result)){result <- 1e12}
    if(!is.numeric(fun_try)){result <- 1e12}
    result
  }
  
  simplexStep <- function(P, FAC) {
    ## Extrapolates by a factor FAC through the face of the simplex across from
    ## the highest (i.e. worst) point.
    worst <- nPOINTS_SIMPLEX
    centr <- apply(P[-worst,,drop=FALSE], 2, mean)
    newpar <- centr*(1-FAC) + P[worst,]*FAC
    newpar
  }
  
  
  if(!file.exists(paste0(state_dir, "/optimization_state.rds"))){
    
    if(!is.null(ini_pop)){
      # predefined population
      cat("Using pre-defined initial population.\n")
      POPULATION <- ini_pop[1:nPOINTS, ]
      spare_samples <- list()
      for(complex in 1:nCOMPLEXES){
        spare_size <- floor((nrow(ini_pop) - nPOINTS) / nCOMPLEXES)
        if(spare_size == 0) {
          spare_samples[[complex]] <- NA
        } else {
          start_id <- nPOINTS + nPOINTS * (complex-1) + 1
          end_id <- start_id + spare_size - 1
          spare_samples[[complex]] <- ini_pop[start_id:end_id, ]
        }
      }
      POP.FITNESS <- numeric(length = nPOINTS)
    } else {
      # sampled population
      # initialize population matrix 
      POPULATION <- matrix(as.numeric(NA), nrow = nPOINTS, ncol = NDIM)
      if (!is.null(names(par)))
        colnames(POPULATION) <- names(par)
      POP.FITNESS <- numeric(length = nPOINTS)
      POPULATION[1,] <- par
      ## generate initial parameter values by random uniform sampling
      finitelower <- ifelse(is.infinite(lower), -(abs(par)+2)*5, lower)
      finiteupper <- ifelse(is.infinite(upper), +(abs(par)+2)*5, upper)
      
      if (control$initsample == "latin") {
        for (i in 1:NDIM) {
          tmp <- seq(finitelower[i], finiteupper[i], length = nPOINTS-1)
          tmp <- jitter(tmp, factor = 2)
          tmp <- pmax(finitelower[i], pmin(finiteupper[i], tmp))
          POPULATION[-1,i] <- sample(tmp)
        }
      } else {
        for (i in 1:NDIM)
          POPULATION[-1,i] <- runif(nPOINTS-1, finitelower[i], finiteupper[i])
      }
    }
    
    ## only store all iterations if requested -- could be big!
    if (!is.finite(MAXIT)) {
      MAXIT <- 10000
      warning("setting maximum iterations to 10000")
    }
    if (returnpop) {
      POP.ALL <- array(as.numeric(NA), dim = c(nPOINTS, NDIM, MAXIT))
      if (!is.null(names(par)))
        dimnames(POP.ALL)[[2]] <- names(par)
    } else {
      POP.ALL <- NA
    }
    POP.FIT.ALL <- matrix(as.numeric(NA), ncol = nPOINTS, nrow = MAXIT)
    BESTMEM.ALL <- matrix(as.numeric(NA), ncol = NDIM, nrow = MAXIT)
    if (!is.null(names(par)))
      colnames(BESTMEM.ALL) <- names(par)
    ## the output object
    obj <- list()
    class(obj) <- c("SCEoptim", class(obj))
    obj$call <- match.call()
    obj$control <- control
    ## initialize timer
    tic <- as.numeric(Sys.time())
    toc <- 0
    ## initialize counters
    funevals <- 0
    
    ## calculate cost for each point in initial population
    # parallelized over multiple instances equal to the number of complexes
    cat("Computing initial population.\n")
    
    # load if previously computed results exist
    if(paste0("parameters_tracker_1.csv") %in% list.files(state_dir)){
      previous_res <- list()
      for(complex in nCOMPLEXES){
        par_res <- read.csv(paste0(state_dir, "/parameters_tracker_", complex, ".csv"))
        loss_res <- read.csv(paste0(state_dir, "/losses_tracker_", complex, ".csv"))
        previous_res[[complex]] <- merge(loss_res[, -c(1:2)], par_res)
      }
      previous_res <- do.call(rbind, previous_res)
      
      # use previous results if already saved
      useable_previous_res <- previous_res[!is.na(previous_res$loss), ]
      useable_previous_res <- useable_previous_res[!duplicated(useable_previous_res[, -c(1:8)]), ]
      # change loss to fit fitness values
      useable_previous_res$loss <- useable_previous_res$loss * fnscale
      idx <- order(useable_previous_res$loss)
      useable_previous_res <- useable_previous_res[idx,]
      useable_previous_res <- useable_previous_res[!is.infinite(useable_previous_res$loss), ]
      
      usable_samples <- min(nrow(useable_previous_res), nPOINTS)
      POPULATION[1:usable_samples, ] <- useable_previous_res[1:usable_samples, -c(1:8)]
      POP.FITNESS[1:usable_samples] <- useable_previous_res[1:usable_samples, "loss"]
      if(usable_samples < nPOINTS) {
        compute_ini_pop <- TRUE
        nPOINTS_COMPLEX <- ceiling((nPOINTS - usable_samples)/nCOMPLEXES)
        start_id <- usable_samples + 1
      } else {
        compute_ini_pop <- FALSE
      }
    } else {
      previous_res <- NULL
      compute_ini_pop <- TRUE
      start_id <- 1
    }
    if(compute_ini_pop){
      init_pop <- function(j){
        k1 <- 1:nPOINTS_COMPLEX
        k2 <- (k1-1) * nCOMPLEXES + j
        k2 <- k2 + start_id - 1
        spare <- 1
        for (i in k2){
          use_old_res <- FALSE
          # check if results already exists
          if(!is.null(previous_res)){
            old_res_avail <- which(!colSums(t(previous_res[,-c(1:8)]) != as.numeric(POPULATION[i, ])))
            if(length(old_res_avail) != 0){
              use_old_res <- TRUE
              if(length(old_res_avail) == 1){
                POP.FITNESS[i] <- previous_res$loss[old_res_avail]  * fnscale
              } else {
                POP.FITNESS[i] <- mean(previous_res$loss[old_res_avail], na.rm=TRUE)  * fnscale
                if(is.nan(POP.FITNESS[i])) POP.FITNESS[i] <- NA
              }
            }
          }
          
          # otherwise compute them
          if(!use_old_res){
            POP.FITNESS[i] <- costFunction(FUN, j, 0, POPULATION[i, ], fnscale)
            
            
            # if NA are produced get different initial sample
            while(POP.FITNESS[i] > 1e11){
              
              # random sample if no spares are available
              if(is.null(ini_pop) | all(is.na(spare_samples[[j]]))){
                POPULATION[i, 1:GenAI_para_dims] <- rnorm(GenAI_para_dims)
                POP.FITNESS[i] <- costFunction(FUN, j, 0, POPULATION[i,], fnscale)
              } 
              # spare sample if available
              if(!all(is.na(spare_samples[[j]]))){
                POPULATION[i, ] <- spare_samples[[j]][spare, ]
                spare <- spare + 1
                POP.FITNESS[i] <- costFunction(FUN, j, 0, POPULATION[i,], fnscale)
                if(nrow(spare_samples[[j]][-c(1:spare), ]) == 0) spare_samples[[j]] <- NA
              }
            }
          }
        }
        return(list(POPULATION, unlist(POP.FITNESS)))
      }
      
      # parallel apply complex steps and merge results
      ini_results <- parallel::mclapply(1:nCOMPLEXES, init_pop)
      cat("**************ini results******************")
      print(str(ini_results))
      for (j in 1:nCOMPLEXES) {
        k1 <- 1:nPOINTS_COMPLEX
        k2 <- (k1-1) * nCOMPLEXES + j
        POPULATION[k2, ] <- ini_results[[j]][[1]][k2, ]
        POP.FITNESS[k2] <- as.numeric(ini_results[[j]][[2]][k2])
      }
    }
    ## sort the population in order of increasing function values
    idx <- order(POP.FITNESS)
    POP.FITNESS <- POP.FITNESS[idx]
    POPULATION <- POPULATION[idx,,drop=FALSE]
    ## store one previous iteration only
    POP.PREV <- POPULATION
    POP.FIT.PREV <- POP.FITNESS
    if (returnpop) {
      POP.ALL[,,1] <- POPULATION
    }
    POP.FIT.ALL[1,] <- POP.FITNESS
    BESTMEM.ALL[1,] <- POPULATION[1,]
    
    ## initiliaize loop variable
    i <- 0
    # save optimizer states
    optimizer_state <- list(
      "tic" = tic,
      "i" = i,
      "POPULATION" = as.matrix(POPULATION),
      "POP.FITNESS" = POP.FITNESS,
      "POP.PREV" = as.matrix(POP.PREV),
      "POP.FIT.PREV" = POP.FIT.PREV,
      "POP.ALL" = as.matrix(POP.ALL),
      "POP.FIT.ALL" = as.matrix(POP.FIT.ALL),
      "BESTMEM.ALL" = as.matrix(BESTMEM.ALL),
      "funevals" = funevals
    )
    saveRDS(optimizer_state, paste0(state_dir, "/optimization_state.rds"))
    
    # store number of iterations
    obj$iterations <- i
    saveRDS(obj, paste0(state_dir, "/last_sce_obj.rds"))
    
    message("Finished computing loss for initial Population")
    return(obj)
  } else {
    # Load states 
    obj <- readRDS(paste0(state_dir, "/last_sce_obj.rds"))
    optimizer_state <- readRDS(paste0(state_dir, "/optimization_state.rds"))
    i <- optimizer_state[["i"]]
    tic <- optimizer_state[["tic"]]
    POPULATION <- optimizer_state[["POPULATION"]]
    POP.FITNESS <- optimizer_state[["POP.FITNESS"]]
    POP.PREV <-  optimizer_state[["POP.PREV"]]
    POP.FIT.PREV <- optimizer_state[["POP.FIT.PREV"]]
    POP.ALL <- optimizer_state[["POP.ALL"]]
    POP.FIT.ALL <- optimizer_state[["POP.FIT.ALL"]]
    BESTMEM.ALL <- optimizer_state[["BESTMEM.ALL"]]
    funevals <- optimizer_state[["funevals"]]
    
    if(any(dim(BESTMEM.ALL) != c(MAXIT, NDIM))){
      BESTMEM.ALL <- unlist(BESTMEM.ALL)
      BESTMEM.ALL <- matrix(BESTMEM.ALL, ncol = NDIM)
    }
    
  }
  
  # set flag and message
  EXITFLAG <- NA
  EXITMSG <- NULL
  ## store best solution from last two iterations
  prevBestVals <- rep(Inf, control$tolsteps)
  prevBestVals[1] <- POP.FITNESS[1]
  
  # optimization loop
  while (i < MAXIT) {
    
    i <- i + 1
    cat("Starting iteration ", i, "/", MAXIT, "\n")
    ## The population matrix POPULATION will now be rearranged into complexes.
    
    ## For each complex ...
    complex_run <- function(j){
      ## construct j-th complex from POPULATION
      
      k1 <- 1:nPOINTS_COMPLEX
      k2 <- (k1-1) * nCOMPLEXES + j
      
      COMPLEX <- POP.PREV[k2,,drop=FALSE]
      COMPLEX_FITNESS <- POP.FIT.PREV[k2]
      
      ## Each complex evolves a number of steps according to the competitive
      ## complex evolution (CCE) algorithm as described in Duan et al. (1992).
      ## Therefore, a number of 'parents' are selected from each complex which
      ## form a simplex. The selection of the parents is done so that the better
      ## points in the complex have a higher probability to be selected as a
      ## parent. The paper of Duan et al. (1992) describes how a trapezoidal
      ## probability distribution can be used for this purpose.
      loop_start_time <- Sys.time()
      for (k in 1:CCEITER) {
        
        # try to skip possible errors without loaing the whole iteration
        check_iter <- try({
          
          ## select simplex by sampling the complex
          
          ## sample points with "trapezoidal" i.e. linear probability
          weights <- rev(ppoints(nPOINTS_COMPLEX))
          ## 'elitism' parameter can give more weight to the better results:
          weights <- weights ^ control$elitism
          LOCATION <- sample(seq(1,nPOINTS_COMPLEX), size = nPOINTS_SIMPLEX,
                             prob = weights)
          
          LOCATION <- sort(LOCATION)
          
          ## construct the simplex
          SIMPLEX <- COMPLEX[LOCATION,,drop=FALSE]
          SIMPLEX_FITNESS <- COMPLEX_FITNESS[LOCATION]
          
          worst <- nPOINTS_SIMPLEX
          
          ## generate new point for simplex
          
          ## first extrapolate by a factor -1 through the face of the simplex
          ## across from the high point,i.e.,reflect the simplex from the high point
          parRef <- simplexStep(SIMPLEX, FAC = -1)
          fitRef <- costFunction(FUN, j, i, parRef, fnscale)
          
          # check if time limit is reached
          keep_running <- TRUE
          if(difftime(Sys.time(), loop_start_time, units = "hours")[[1]] > time_limit){
            keep_running <- FALSE
          }
          
          ## check the result
          if (fitRef <= SIMPLEX_FITNESS[1] & keep_running) {
            ## gives a result better than the best point,so try an additional
            ## extrapolation by a factor 2
            parRefEx <- simplexStep(SIMPLEX, FAC = -2)
            cat("Iteration run time:", difftime(Sys.time(), loop_start_time, units = "hours")[[1]], "h")
            fitRefEx <- costFunction(FUN, j, i, parRefEx, fnscale)
            
            if (fitRefEx < fitRef) {
              SIMPLEX[worst,] <- parRefEx
              SIMPLEX_FITNESS[worst] <- fitRefEx
              ALGOSTEP <- 'reflection and expansion'
            } else {
              SIMPLEX[worst,] <- parRef
              SIMPLEX_FITNESS[worst] <- fitRef
              ALGOSTEP <- 'reflection'
            }
          } else if (fitRef >= SIMPLEX_FITNESS[worst-1] & keep_running) {
            ## the reflected point is worse than the second-highest, so look
            ## for an intermediate lower point, i.e., do a one-dimensional
            ## contraction
            parCon <- simplexStep(SIMPLEX, FAC = -0.5)
            fitCon <- costFunction(FUN, j, i, parCon, fnscale)
            
            if (fitCon < SIMPLEX_FITNESS[worst]) {
              SIMPLEX[worst,] <- parCon
              SIMPLEX_FITNESS[worst] <- fitCon
              ALGOSTEP <- 'one dimensional contraction'
            } else {
              ## can't seem to get rid of that high point, so better contract
              ## around the lowest (best) point
              SIMPLEX <- (SIMPLEX + rep(SIMPLEX[1,], each=nPOINTS_SIMPLEX)) / 2
              for (k in 2:nrow(SIMPLEX)){
                SIMPLEX_FITNESS[k] <- costFunction(FUN, j, i, SIMPLEX[k,], fnscale)
              }
              ALGOSTEP <- 'multiple contraction'
            }
          } else if (keep_running) {
            ## if better than second-highest point, use this point
            SIMPLEX[worst,] <- parRef
            SIMPLEX_FITNESS[worst] <- fitRef
            ALGOSTEP <- 'reflection'
          }
          
          message(ALGOSTEP)
          
          ## replace the simplex into the complex
          COMPLEX[LOCATION,] <- SIMPLEX
          COMPLEX_FITNESS[LOCATION] <- SIMPLEX_FITNESS
          
          ## sort the complex
          idx <- order(COMPLEX_FITNESS)
          COMPLEX_FITNESS <- COMPLEX_FITNESS[idx]
          COMPLEX <- COMPLEX[idx,,drop=FALSE]
        })
        
        
        # check if time limit is reached
        if(difftime(Sys.time(), loop_start_time, units = "hours")[[1]] > time_limit){
          break
        } else {
          keep_running <- TRUE
        }
        
      }
      return(list(k2, COMPLEX, unlist(COMPLEX_FITNESS)))
    }
    # parallel apply complex steps and merge results
    cat("**************complex results******************")
    complex_results <- parallel::mclapply(1:nCOMPLEXES, complex_run)
    cat("\n parallel compex results structure:\n")
    print(str(complex_results))
    for (j in 1:nCOMPLEXES) {
      # ## replace the complex back into the population
      k2 <- complex_results[[j]][[1]]
      COMPLEX <- complex_results[[j]][[2]]
      COMPLEX_FITNESS <- complex_results[[j]][[3]]
      POPULATION[k2,] <- COMPLEX
      POP.FITNESS[k2] <- COMPLEX_FITNESS
    }
    
    ## At this point, the population was divided in several complexes, each of which
    ## underwent a number of iteration of the simplex (Metropolis) algorithm. Now,
    ## the points in the population are sorted, the termination criteria are checked
    ## and output is given on the screen if requested.
    
    ## sort the population
    idx <- order(POP.FITNESS)
    POP.FITNESS <- POP.FITNESS[idx]
    POPULATION <- POPULATION[idx,,drop=FALSE]
    if (returnpop) {
      POP.ALL[,,i] <- POPULATION
    }
    POP.FIT.ALL[i,] <- POP.FITNESS
    BESTMEM.ALL[i,] <- POPULATION[1,]
    
    curBest <- POP.FITNESS[1]
    
    ## end the optimization if one of the stopping criteria is met
    
    prevBestVals <- c(curBest, head(prevBestVals, -1))
    reltol <- control$reltol
    if (all(abs(diff(prevBestVals)) <= reltol * (abs(curBest)+reltol))) {
      EXITMSG <- 'Change in solution over [tolsteps] less than specified tolerance (reltol).'
      EXITFLAG <- 0
    }
    
    ## give user feedback on screen if requested
    if (trace >= 1) {
      if (i == 1) {
        message(' Nr Iter  Nr Fun Eval    Current best function    Current worst function')
      }
      if ((i %% control$REPORT == 0) || (!is.na(EXITFLAG)))
      {
        message(sprintf(' %5.0f     %5.0f             %12.6g              %12.6g',
                        i, funevals, min(POP.FITNESS), max(POP.FITNESS)))
        if (trace >= 2)
          message("parameters: ", toString(signif(POPULATION[1,], 3)))
      }
    }
    
    if (!is.na(EXITFLAG))
      break
    
    if ((i >= control$maxit) || (funevals >= control$maxeval)) {
      EXITMSG <- 'Maximum number of function evaluations or iterations reached.'
      EXITFLAG <- 1
      break
    }
    
    toc <- as.numeric(Sys.time()) - tic
    if (toc > control$maxtime) {
      EXITMSG <- 'Exceeded maximum time.'
      EXITFLAG <- 2
      break
    }
    
    ## go to next iteration
    POP.PREV <- POPULATION
    POP.FIT.PREV <- POP.FITNESS
    
    # save optimizer states
    optimizer_state <- list(
      "tic" = tic,
      "i" = i,
      "POPULATION" = POPULATION,
      "POP.FITNESS" = POP.FITNESS,
      "POP.PREV" = POP.PREV,
      "POP.FIT.PREV" = POP.FIT.PREV,
      "POP.ALL" = POP.ALL,
      "POP.FIT.ALL" = POP.FIT.ALL,
      "BESTMEM.ALL" = BESTMEM.ALL,
      "funevals" = funevals
    )
    ## return solution
    obj$par <- POPULATION[1,]
    obj$value <- POP.FITNESS[1]
    obj$convergence <- EXITFLAG
    obj$message <- EXITMSG
    
    ## store number of function evaluations
    obj$counts <- funevals
    ## store number of iterations
    obj$iterations <- i
    ## store the amount of time taken
    obj$time <- toc
    
    saveRDS(optimizer_state, paste0(state_dir, "/optimization_state.rds"))
    saveRDS(obj, paste0(state_dir, "/last_sce_obj.rds"))
    
    # save also in iteration specific order
    iter_result_dir <- paste0(state_dir, "/iteration_results")
    dir.create(iter_result_dir, showWarnings = FALSE)
    saveRDS(optimizer_state, paste0(iter_result_dir, "/iter_", i, " _optimization_state.rds"))
    
    if(restart_after_n_iter == 1){
      restart <- TRUE
    } else if(i %% restart_after_n_iter == 0){
      restart <- TRUE
    } else {
      restart <- FALSE
    }
    if(restart) if(i < MAXIT) return(obj)
  }
  
  if (trace >= 1)
    message(EXITMSG)
  
  ## return solution
  obj$par <- POPULATION[1,]
  obj$value <- POP.FITNESS[1]
  obj$convergence <- EXITFLAG
  obj$message <- EXITMSG
  
  ## store number of function evaluations
  obj$counts <- funevals
  ## store number of iterations
  obj$iterations <- i
  ## store the amount of time taken
  obj$time <- toc
  
  if (returnpop) {
    ## store information on the population at each iteration
    obj$POP.ALL <- POP.ALL[,,1:i]
    dimnames(obj$POP.ALL)[[3]] <- paste("iteration", 1:i)
  }
  obj$POP.FIT.ALL <- POP.FIT.ALL[1:i,]
  obj$BESTMEM.ALL <- BESTMEM.ALL[1:i,]
  saveRDS(obj, paste0(state_dir, "/last_sce_obj.rds"))
  
  return(obj)
}


# nml utilities --------------------------------------------------------------------------


write_mhm_nml  <-	function(glm_nml,file){
  sink(file)
  cat("! Emacs: -*- mode: f90 -*-\n",
      "!global_parameters\n",
      "!PARAMETER   lower_bound  upper_bound  value  FLAG  SCALING\n",
      "!interception)\n")
  print(glm_nml)
  sink()
}

read_nml  <-	function(nml_file){
  # skip all commented lines, return all variables and associated values
  # requires NO return line variables (all variables must be completely defined on a single line)
  c <- file(nml_file,"r")
  fileLines <- readLines(c)
  close(c)
  lineStart	<-	substr(fileLines,1,1)
  # ignore comment lines or empty lines
  ignoreLn	<-	lineStart=='!' | fileLines==""
  lineStart	<-	lineStart[!ignoreLn]
  fileLines	<-	fileLines[!ignoreLn]
  # find all lines which start with "&" * requires FIRST char to be value
  
  lineIdx		<- seq(1,length(lineStart))
  blckOpen	<-	lineIdx[lineStart=="&"]
  blckClse	<-	lineIdx[lineStart=="/"]
  
  nml <- list()
  for (i in seq_len(length(blckOpen))){
    blckName   <-	substr(fileLines[blckOpen[i]],
                         2, nchar(fileLines[blckOpen[i]]))
    blckName   <- gsub("\\s", "", blckName)
    oldNms	   <-	names(nml)
    nml[[i]]   <-	list()
    names(nml) <-	c(oldNms,blckName)
    
    carryover <- ''
    
    for (j in (blckOpen[i]+1):(blckClse[i]-1)){
      
      textLine	<-	paste(carryover,
                        gsub("\t", "", gsub(" ", "", fileLines[j])), sep = '')
      
      if(substr(textLine, 1, 1) != '!'){
        # Add a check here, sometimes, if there is a hanging comma,
        #and only sometimes that means add next row
        if(substr(textLine, nchar(textLine), nchar(textLine)) == ',' &&
           j+1 <= length(fileLines) &&
           !any(grep("=", fileLines[j + 1])) &&
           !any(grep("/", fileLines[j + 1]))){
          
          carryover = textLine
          next
        }else{
          carryover = ''
        }
        # else, line is commented out
        lineVal	  <-	buildVal(textLine, lineNum = j, blckName)
        nml[[i]]	<-	c(nml[[i]], lineVal)
      }
    }
  }
  nml <- .nml(nml)
  return(nml)
}

print.nml <- function(x, ...){
  glm_nml <- x
  for (i in seq_len(length(names(glm_nml)))){ # these are the blocks
    blckNm  <-	names(glm_nml)[i]
    cat("&")
    cat(blckNm)
    cat('\n')
    blckList	<-	glm_nml[[i]]
    for (j in seq_len(length(names(blckList)))){
      cat('   ')
      cat(names(blckList)[j])
      cat(' = ')
      if (length(blckList[[j]])>1){
        if (is.logical(blckList[[j]])){
          charText <- to.glm_boolean(blckList[[j]])
        } else {
          charText <- c(blckList[[j]])
        }
        writer	<-	paste(charText,collapse=', ')
      } else if (is.character(blckList[[j]])) {
        charText <- strsplit(blckList[[j]],',')
        writer <- paste(c("'",paste(c(charText[[1]]),collapse="','"),"'"),collapse='')
      } else if (is.logical(blckList[[j]])){
        writer <- to.glm_boolean(blckList[[j]])
      } else {
        writer <- blckList[[j]]
      }
      cat(writer)
      cat('\n')
    }
    cat('/\n')
  }
}

buildVal	<-	function(textLine, lineNum, blckName){
  #-----function appends nml list with new values-----
  # remove all text after comment string
  textLine	<-	strsplit(textLine,'!')[[1]][1]
  
  if (!any(grep("=", textLine))){
    stop(c("no hanging lines allowed in .nml, used ",textLine,'.\nSee line number:',lineNum,' in "&',blckName,'" section.'))
  }
  params	<-	strsplit(textLine,"=") # break text at "="
  parNm	  <-	params[[1]][1]
  parVl	  <-	params[[1]][2]
  # figure out what parval is...if string, remove quotes and keep as string
  # ***for boolean text, use "indentical" so that 0!= FALSE
  # can be: string, number, comma-sep-numbers, or boolean
  
  # special case for date:
  if (is.na(parVl)){
    stop('Empty values after "', textLine, '" on line ', lineNum,
         '. \nPerhaps the values are on the next line?', call. = FALSE)
  }
  if (nchar(parVl>17) & substr(parVl,14,14)==':' & substr(parVl,17,17)==':'){
    parVl<-paste(c(substr(parVl,1,11),' ',substr(parVl,12,nchar(parVl))),collapse='')
  }
  if (any(grep("'",parVl))){
    
    parVl	<-	gsub("'","",parVl)
  }else if (any(grep("\"",parVl))){
    parVl  <-	gsub("\"","",parVl)
  }else if (isTRUE(grepl(".true.",parVl) || grepl(".false.",parVl))){
    logicals <- unlist(strsplit(parVl,","))
    parVl <- from.glm_boolean(logicals)
  }else if (any(grep(",",parVl))){	# comma-sep-nums
    parVl	<-	c(as.numeric(unlist(strsplit(parVl,","))))
  }else {	# test for number
    parVl	<-	as.numeric(parVl)
  }
  lineVal	<-	list(parVl)
  names(lineVal)	<-	parNm
  return(lineVal)
}

from.glm_boolean <- function(values){
  logicals <- sapply(values, FUN = function(x){
    if (!isTRUE(grepl(".true.", x) || grepl(".false.", x))){
      stop(x, ' is not a .true. or .false.; conversion to TRUE or FALSE failed.',
           call. = FALSE)
    }
    return(ifelse(isTRUE(grepl(".true.", x)), TRUE, FALSE))
  })
  return(as.logical(logicals))
}

to.glm_boolean <- function(values){
  val.logical <- values
  values[val.logical] <- '.true.'
  values[!val.logical] <- '.false.'
  return(values)
}

findBlck	<-	function(nml,argName){
  
  # test for argName being a string
  if (!is.character(argName)){stop(c("parameter name must be a string"))}
  fau <- " "
  fault.string <- rep(fau,1000) # names fault matrix, only returned when empty match
  blockNames	<-	names(nml)
  blckI	<-	c()
  for (i in seq_len(length(blockNames))){
    if (any(argName %in% names(nml[[i]]))){
      blckI	<- c(blckI,i)
    } else {
      one.i <- which(fault.string==fau)[1]
      fault.string[one.i:(one.i+length(names(nml[[i]]))-1)]=names(nml[[i]])
    }
    
  }
  fault.string <- fault.string[!fault.string==fau] # is empty if found
  # test to see if a block match was made
  if (is.null(blckI)){stop(c("parameter name ",argName," not found in nml. Possible names:",paste(fault.string,collapse=', ')))}
  return(blckI)
}

setnmlList <- function(glm_nml,arg_list){
  if (!is.list(arg_list)){stop("arg_list must be a list")}
  if (any(nchar(names(arg_list)) == 0) | length(names(arg_list)) == 0){
    stop('arg_list must be a named list')
  }
  arg_names  <-	names(arg_list)
  for (i in seq_len(length(arg_names))){
    glm_nml <- set_nml(glm_nml,arg_name=arg_names[i],arg_val=arg_list[[i]])
  }
  return(glm_nml)
}


is_nml_file <- function(nml_file){
  is_nml <- FALSE
  fl_ext <- tail(strsplit(nml_file, "\\.")[[1]],1)
  if (fl_ext == 'nml'){
    is_nml <- TRUE
  }
  return(is_nml)
}

what_ascii <- function(file){
  response <- capture.output(showNonASCIIfile(file))
  return(response)
}

ascii_only <- function(file){
  response <- what_ascii(file)
  if (length(response) > 0){
    return(FALSE)
  } else {
    return(TRUE)
  }
}


get_block <- function(glm_nml, arg_name, warn=TRUE){
  arg_split = strsplit(arg_name,'::')[[1]]
  if (length(arg_split) > 1){
    blck = arg_split[1]
    arg_name = get_arg_name(arg_name)
  } else{
    blck	<-	findBlck(glm_nml,arg_name)
  }
  if(length(blck) > 1){
    if(warn){
      warning(arg_name, " found in ", 
              paste(names(glm_nml[blck]), collapse=' & '), 
              ", returning the first. Try ",
              names(glm_nml[blck])[1],"::",arg_name, " for explicit match")
    }
    blck = blck[1]
  }
  return(blck)
}

get_arg_name <- function(arg_name){
  arg_split = strsplit(arg_name,'::')[[1]]
  
  if (length(arg_split) > 1){
    blck = arg_split[1]
    arg_name = arg_split[2]
  }
  return(arg_name)
}

.nml <- function(list_obj){
  nml <- list_obj
  class(nml) <- "nml"
  invisible(nml)
}


# Quality Criteria -----------------------------------------------------------------------
# Functions taken from https://github.com/hzambran/hydroGOF


NSE.default <- function (sim, obs, na.rm=TRUE, fun=NULL, ..., 
                         epsilon.type=c("none", "Pushpalatha2012", "otherFactor", "otherValue"), 
                         epsilon.value=NA){ 
  
  if ( is.na(match(class(sim), c("integer", "numeric", "ts", "zoo", "xts"))) |
       is.na(match(class(obs), c("integer", "numeric", "ts", "zoo", "xts")))
  ) stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo', 'xts')")      
  
  epsilon.type <- match.arg(epsilon.type)  
  
  # index of those elements that are present both in 'sim' and 'obs' (NON- NA values)
  vi <- valindex(sim, obs)
  
  if (length(vi) > 0) {	 
    # Filtering 'obs' and 'sim', selecting only those pairs of elements 
    # that are present both in 'x' and 'y' (NON- NA values)
    obs <- obs[vi]
    sim <- sim[vi]
    
    if (!is.null(fun)) {
      fun1 <- match.fun(fun)
      new  <- preproc(sim=sim, obs=obs, fun=fun1, ..., 
                      epsilon.type=epsilon.type, epsilon.value=epsilon.value)
      sim  <- new[["sim"]]
      obs  <- new[["obs"]]
    } # IF end     
    
    denominator <- sum( (obs - mean(obs))^2 )
    
    if (denominator != 0) {      
      NS <- 1 - ( sum( (obs - sim)^2 ) / denominator )     
    } else {
      NS <- NA
      warning("'sum((obs - mean(obs))^2)=0' => it is not possible to compute 'NSE'")  
    } 
  } else {
    NS <- NA
    warning("There are no pairs of 'sim' and 'obs' without missing values !")
  } # ELSE end
  
  return(NS)
  
}

NSE.matrix <- function(sim, obs, na.rm=TRUE, fun=NULL, ..., 
                       epsilon.type=c("none", "Pushpalatha2012", "otherFactor", "otherValue"), 
                       epsilon.value=NA){ 
  
  # Checking that 'sim' and 'obs' have the same dimensions
  if ( all.equal(dim(sim), dim(obs)) != TRUE )
    stop( paste("Invalid argument: dim(sim) != dim(obs) ( [", 
                paste(dim(sim), collapse=" "), "] != [", 
                paste(dim(obs), collapse=" "), "] )", sep="") )
  
  NS <- rep(NA, ncol(obs))       
  
  NS <- sapply(1:ncol(obs), function(i,x,y) { 
    NS[i] <- NSE.default( x[,i], y[,i], na.rm=na.rm, fun=fun, ..., 
                          epsilon.type=epsilon.type, epsilon.value=epsilon.value)
  }, x=sim, y=obs )    
  
  names(NS) <- colnames(obs)
  
  return(NS)
  
}

NSE.data.frame <- function(sim, obs, na.rm=TRUE, fun=NULL, ..., 
                           epsilon.type=c("none", "Pushpalatha2012", "otherFactor", "otherValue"), 
                           epsilon.value=NA){ 
  
  sim <- as.matrix(sim)
  obs <- as.matrix(obs)
  
  NSE.matrix(sim, obs, na.rm=na.rm, fun=fun, ..., 
             epsilon.type=epsilon.type, epsilon.value=epsilon.value)
}
NSE <-function(sim, obs, ...) UseMethod("NSE")

KGE.default <- function(sim, obs, s=c(1,1,1), na.rm=TRUE, 
                        method=c("2009", "2012"), out.type=c("single", "full"), 
                        fun=NULL, ...,
                        epsilon.type=c("none", "Pushpalatha2012", "otherFactor", "otherValue"), 
                        epsilon.value=NA) { 
  
  # If the user provided a value for 's'
  if (!identical(s, c(1,1,1)) )  {
    if ( length(s) != 3 ) stop("Invalid argument: lenght(s) must be equal to 3 !")
    if ( sum(s) != 1 )    stop("Invalid argument: sum(s) must be equal to 1.0 !")
  } # IF end
  
  method   <- match.arg(method)
  out.type <- match.arg(out.type)  
  
  if ( is.na(match(class(sim), c("integer", "numeric", "ts", "zoo"))) |
       is.na(match(class(obs), c("integer", "numeric", "ts", "zoo")))
  ) stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")      
  
  vi <- valindex(sim, obs)
  
  if (length(vi) > 0) {
    
    obs <- as.numeric(obs[vi])
    sim <- as.numeric(sim[vi])
    
    if (!is.null(fun)) {
      fun1 <- match.fun(fun)
      new  <- preproc(sim=sim, obs=obs, fun=fun1, ..., 
                      epsilon.type=epsilon.type, epsilon.value=epsilon.value)
      sim  <- new[["sim"]]
      obs  <- new[["obs"]]
    } # IF end
    
    # Mean values
    mean.sim <- mean(sim, na.rm=na.rm)
    mean.obs <- mean(obs, na.rm=na.rm)
    
    # Standard deviations
    sigma.sim <- sd(sim, na.rm=na.rm)
    sigma.obs <- sd(obs, na.rm=na.rm)
    
    # Pearson product-moment correlation coefficient
    r     <- rPearson(sim, obs)
    
    # Alpha is a measure of relative variability between simulated and observed values (See Ref1)
    Alpha <- sigma.sim / sigma.obs
    
    # Beta is the ratio between the mean of the simulated values to the mean of observations
    Beta <- mean.sim / mean.obs
    
    # CV.sim is the coefficient of variation of the simulated values [dimensionless]
    # CV.obs is the coefficient of variation of the observations [dimensionless]
    CV.sim <- sigma.sim / mean.sim
    CV.obs <- sigma.obs / mean.obs
    
    # Gamma is the variability ratio, which is used instead of Alpha (See Ref2)
    Gamma <- CV.sim / CV.obs
    
    # Variability ratio depending on 'method'
    if(method=="2012") {
      vr     <- Gamma
      vr.stg <- "Gamma"
    } else {
      vr     <- Alpha
      vr.stg <- "Alpha"
    } # ELSE end
    
    # KGE Computation
    if ( (mean.obs != 0) | (sigma.obs != 0) ) {
      KGE <- 1 - sqrt( (s[1]*(r-1))^2 + (s[2]*(vr-1))^2 + (s[3]*(Beta-1))^2 )
    } else {
      if ( mean.obs != 0)  warning("Warning: 'mean(obs)==0'. Beta = Inf")
      if ( sigma.obs != 0) warning("Warning: 'sd(obs)==0'. ", vr.stg, " = Inf")
      KGE <- NA
    } # ELSE end  
    
  } else {
    r    <- NA
    Beta <- NA
    vr   <- NA
    if(method=="2012") {
      vr.stg <- "Gamma"
    } else vr.stg <- "Alpha" 
    KGE <- NA
    warning("There are no pairs of 'sim' and 'obs' without missing values !")
  } # ELSE end
  
  if (out.type=="single") {
    out <- KGE
  } else {
    out <- list(KGE.value=KGE, KGE.elements=c(r, Beta, vr))
    names(out[[2]]) <- c("r", "Beta", vr.stg)
  } # ELSE end    
  
  return(out)
  
}


KGE.matrix <- function (sim, obs, s=c(1,1,1), na.rm=TRUE, 
                        method=c("2009", "2012"), out.type=c("single", "full"), 
                        fun=NULL, ...,
                        epsilon.type=c("none", "Pushpalatha2012", "otherFactor", "otherValue"), 
                        epsilon.value=NA) { 
  
  # Checking that 'sim' and 'obs' have the same dimensions
  if ( all.equal(dim(sim), dim(obs)) != TRUE )
    stop( paste("Invalid argument: dim(sim) != dim(obs) ( [", 
                paste(dim(sim), collapse=" "), "] != [", 
                paste(dim(obs), collapse=" "), "] )", sep="") )
  
  # If the user provided a value for 's'
  if (!all.equal(s, c(1,1,1)) )  {
    if ( length(s) != 3 ) stop("Invalid argument: lenght(s) must be equal to 3 !")
    if ( sum(s) != 1 )    stop("Invalid argument: sum(s) must be equal to 1.0 !")
  } # IF end
  
  method   <- match.arg(method)
  out.type <- match.arg(out.type) 
  
  ifelse(method=="2012", vr.stg <- "Gamma", vr.stg <- "Alpha")
  
  KGE                <- rep(NA, ncol(obs))       
  elements           <- matrix(NA, nrow=3, ncol=ncol(obs))
  rownames(elements) <- c("r", "Beta", vr.stg)
  colnames(elements) <- colnames(obs)
  
  if (out.type=="single") {
    out <- sapply(1:ncol(obs), function(i,x,y) { 
      KGE[i] <- KGE.default( x[,i], y[,i], s=s, na.rm=na.rm, 
                             method=method, out.type=out.type, 
                             fun=fun, ..., epsilon.type=epsilon.type, 
                             epsilon.value=epsilon.value )
    }, x=sim, y=obs )  
    names(out) <- colnames(obs) 
  } else { out <- lapply(1:ncol(obs), function(i,x,y) { 
    KGE.default( x[,i], y[,i], s=s, na.rm=na.rm, method=method, 
                 out.type=out.type, fun=fun, ..., 
                 epsilon.type=epsilon.type, 
                 epsilon.value=epsilon.value )
  }, x=sim, y=obs ) 
  for (i in 1:length(out) ) {
    KGE[i] <- out[[i]][[1]]
    elements[,i] <- as.numeric(out[[i]][[2]])
  } # FOR end 
  out <- list(KGE.value=KGE, KGE.elements=elements)
  } # ELSE end                     
  
  return(out)
  
}

KGE.data.frame <- function (sim, obs, s=c(1,1,1), na.rm=TRUE, 
                            method=c("2009", "2012"), out.type=c("single", "full"), 
                            fun=NULL, ...,
                            epsilon.type=c("none", "Pushpalatha2012", "otherFactor", "otherValue"), 
                            epsilon.value=NA) { 
  
  sim <- as.matrix(sim)
  obs <- as.matrix(obs)
  
  method   <- match.arg(method)
  out.type <- match.arg(out.type) 
  
  KGE.matrix(sim, obs, s=s, na.rm=na.rm, method=method, out.type=out.type, 
             fun=fun, ..., epsilon.type=epsilon.type, epsilon.value=epsilon.value)
  
}
KGE <- function(sim, obs, ...) UseMethod("KGE")

valindex.default <- function(sim, obs, ...) {  
  
  if ( length(obs) != length(sim) ) {
    stop( "Invalid argument: 'length(sim) != length(obs)' !! (", length(sim), "!=", length(obs), ") !!" )
  } else { 
    index <- which(!is.na(sim) & !is.na(obs))
    if (length(index)==0) warning("'sim' and 'obs' are empty or they do not have any common pair of elements with data !!")
    return( index  )
  } # ELSE end
  
} 

valindex.matrix <- function(sim, obs, ...) { 
  
  # Checking that 'sim' and 'obs' have the same dimensions
  if ( all.equal(dim(sim), dim(obs)) != TRUE ) {
    stop( paste("Invalid argument: dim(sim) != dim(obs) ( [", 
                paste(dim(sim), collapse=" "), "] != [", 
                paste(dim(obs), collapse=" "), "] )", sep="") )
  } else  
    return ( !is.na( sim) & !is.na(obs) )
  
}
valindex <- function(sim, obs, ...) UseMethod("valindex")


rPearson.default <- function(sim, obs, fun=NULL, ..., 
                             epsilon.type=c("none", "Pushpalatha2012", "otherFactor", "otherValue"), 
                             epsilon.value=NA) {
  
  if ( is.na(match(class(sim), c("integer", "numeric", "ts", "zoo"))) |
       is.na(match(class(obs), c("integer", "numeric", "ts", "zoo")))
  ) stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")
  
  vi <- valindex(sim, obs)
  
  if (length(vi) > 0) {
    
    obs <- as.numeric(obs[vi])
    sim <- as.numeric(sim[vi])
    
    if (!is.null(fun)) {
      fun1 <- match.fun(fun)
      new  <- preproc(sim=sim, obs=obs, fun=fun1, ..., 
                      epsilon.type=epsilon.type, epsilon.value=epsilon.value)
      sim  <- new[["sim"]]
      obs  <- new[["obs"]]
    } # IF end
    
    rPearson <- cor(sim, obs, method="pearson", use="pairwise.complete.obs")      
    # if 'sim' and 'obs' were matrixs or data.frame, then the correlation
    # between observed and simulated values for each variable is given by the diagonal of 'r.Pearson' 
    
    #if ( is.matrix(r.Pearson) | is.data.frame(r.Pearson) ) {
    #r.Pearson        <- diag(r.Pearson)
    #}
    
  } else {
    rPearson <- NA
    warning("There are no pairs of 'sim' and 'obs' without missing values !")
  } # ELSE end
  
  return(rPearson)
  
} 

rPearson.matrix <- function(sim, obs, na.rm=TRUE, fun=NULL, ..., 
                            epsilon.type=c("none", "Pushpalatha2012", "otherFactor", "otherValue"), 
                            epsilon.value=NA){
  
  rPearson <- rep(NA, ncol(obs))       
  
  rPearson <- sapply(1:ncol(obs), function(i,x,y) { 
    rPearson[i] <- rPearson.default( x[,i], y[,i], na.rm=na.rm, fun=fun, ..., 
                                     epsilon.type=epsilon.type, epsilon.value=epsilon.value)
  }, x=sim, y=obs )            
  
  return(rPearson)
  
} 

rPearson.data.frame <- function(sim, obs, na.rm=TRUE, fun=NULL, ..., 
                                epsilon.type=c("none", "Pushpalatha2012", "otherFactor", "otherValue"), 
                                epsilon.value=NA){
  
  sim <- as.matrix(sim)
  obs <- as.matrix(obs)
  
  rPearson.matrix(sim, obs, na.rm=na.rm, fun=fun, ..., 
                  epsilon.type=epsilon.type, epsilon.value=epsilon.value)        
  
} 
rPearson <-function(sim, obs, ...) UseMethod("rPearson")
