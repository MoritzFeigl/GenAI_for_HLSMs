# Reset run with updated initial population


grid <- data.frame("exp" = c(0, 1, 2),
                   "run" = c(1, 1, 1))

for(i in 1:nrow(grid)){
  # define run name to create environment automatically
  
  experiment <- grid[i, "exp"]
  run <- grid[i, "run"]
  run_name <- paste0("exp", experiment, "-run", run)
  main_path <- "/gpfs/data/fs71468/GenAI_para_runs/"
  
  # Environment setup ----------------------------------------------------------------------
  # create run and result folders
  run_dir <- paste0(main_path, "runs/", run_name)
  result_dir <- paste0(main_path, "results/", run_name)
  code_dir <- paste0(main_path, "code/")
  data_dir <- paste0(main_path, "data/")
  setwd(run_dir)
  
  # load utils
  source(paste0(code_dir, "02_utility_functions/01_GenAI_para_utility_functions.R"))
  source(paste0(code_dir, "02_utility_functions/02_vae_generators.R"))
  
  # define experiment configuration
  source(paste0(code_dir, 
                "06_run_experiments/experiment_configs/experiment_", experiment, "_setup.R"))
  
  # parallel setups is bound to number of complexes
  n_parallel_setups <- ncomplexes
  
  # Get all previous results and combine to new population
  para1 <- read.csv(paste0(result_dir, "/parameters_tracker_1.csv"))
  para2 <- read.csv(paste0(result_dir, "/parameters_tracker_2.csv"))
  paras <- rbind(para1, para2)
  
  loss1 <- read.csv(paste0(result_dir, "/losses_tracker_1.csv"))
  loss2 <- read.csv(paste0(result_dir, "/losses_tracker_2.csv"))
  losses <- rbind(loss1, loss2)
  
  new_pop_base <- merge(losses[, c("complex", "time", "loss")], paras[, -1], by=c("complex", "time"))
  
  # remove NA losses
  new_pop_base <- new_pop_base[!is.na(new_pop_base$loss), ]
  # remove duplicates
  new_pop_base <- new_pop_base[-which(duplicated(new_pop_base[, -c(1:3)])), ]
  # order by loss
  new_pop_base <- new_pop_base[order(new_pop_base$loss, decreasing = TRUE), ]
  
  # Define new Population
  new_pop <- new_pop_base[, -c(1:3)]
  write.csv(new_pop, paste0(run_dir, "/initial_population.csv"), row.names = FALSE)
  
  # Rename old results and population file
  file.rename(paste0(run_dir, "/initial_population.csv"), paste0(run_dir, "/old_inipop_initial_population.csv"))
  result_files <- list.files(result_dir)
  new_names <- paste0("old_inipop_", result_files)
  file.rename(paste0(result_dir, "/", result_files), paste0(result_dir, "/", new_names))
  
  # restart run
  bash_file <- paste0("/gpfs/data/fs71468/GenAI_para_runs/runs/exp",
                      experiment, "-run", run,
                      "/01_experiment_start.sh")
  system(paste0("sbatch ", bash_file))
  
}