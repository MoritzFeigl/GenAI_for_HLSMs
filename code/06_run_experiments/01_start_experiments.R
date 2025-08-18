# GenAI_para-ET-mHM project
# start optimization runs

# Create run scripts for GenAI_para-mHM experiments
main_path <- "/gpfs/data/fs71468/GenAI_para_runs/"
nruns <- 5
experiments <- c(0, 1, 2)

# create main directories
code_dir <- paste0(main_path, "code")
dir.create(paste0(main_path, "runs"), showWarnings = FALSE)
dir.create(paste0(main_path, "results"), showWarnings = FALSE)

for(experiment in experiments){
  for(run in 1:nruns){
    # define run name to create environment automatically
    run_name <- paste0("exp", experiment, "-run", run)
    
    # main directory names
    run_dir <- paste0(main_path, "runs/", run_name)
    
    # create dirs
    dir.create(run_dir, showWarnings = FALSE)

    # adapt GenAI_para run script
    u03 <- readLines(paste0(code_dir, "/02_utility_functions/03_mhm_GenAI_para_runs.R"))
    u03[grep("experiment <-", u03)] <- paste0("experiment <- ", experiment)
    u03[grep("run_name <-", u03)] <- paste0('run_name <- "', run_name, '"')
    u03[grep("main_path <-", u03)] <- paste0('main_path <- "', main_path, '"')
    writeLines(u03, paste0(run_dir, "/02_mhm_GenAI_para_runs.R"))
    
    # get sbatch script and save as experiment start bash
    u06 <- readLines(paste0(code_dir, "/02_utility_functions/06_start_GenAI_para_mhm.sh"))
    u06[grep("#SBATCH -J ", u06)] <- paste0("#SBATCH -J ", run_name)
    u06[grep("cd ", u06)] <- paste0("cd ", run_dir)
    writeLines(u06, paste0(run_dir, "/01_experiment_start.sh"))
    
    # start optimization
    setwd(run_dir)
    system("chmod 775 01_experiment_start.sh")
    system("sbatch 01_experiment_start.sh")
  }
}
