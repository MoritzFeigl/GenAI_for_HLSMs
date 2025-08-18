
keep_running <- TRUE
start_script <- Sys.time()
while(keep_running){
  
  # Get information on current running jobs
  running <- system("squeue --nohead --format %i-%j-%t -u mfeigl", intern = TRUE)
  FormatRunInfo <- function(squeue_output){
    info <- unlist(strsplit(squeue_output, "-"))
    data.frame("id" = info[1],
               "experiment" = gsub("exp", "", info[2]),
               "run" = gsub("run", "", info[3]),
               "status" = info[4])
  }
  runs <- do.call(rbind, lapply(running, FormatRunInfo))
  runs <- runs[order(runs$experiment),]
  runs$curr_iter <- runs$it_sinc_last_val_2 <- runs$it_sinc_last_val_1 <- runs$SPAEF <- runs$KGE <- 
    runs$lNSE  <- runs$NSE  <- runs$loss <- runs$time  <- runs$iteration <- NA
  
  # merge with jobs that should be running
  all_jobs <- data.frame("experiment" = rep(c(0, 1, 2), each = 5),
                         "run" = rep(1:5, times = 3))
  runs <- merge(all_jobs, runs, by = c("experiment", "run"), all = TRUE)
  
  for(run in 1:nrow(runs)){
    
    restart <- FALSE
    
    if(is.na(runs$id[run])) restart <- TRUE
    
    # check if CG and restart in necessary
    try({if(runs[run, "status"] == "CG") restart <- TRUE}, silent = TRUE)
    try({if(runs[run, "status"] == "PD") next}, silent = TRUE)
    
    # check if run is frozen or otherwisely stopped
      try({
        # get loss and NSE
        result_1 <- read.csv(paste0("/gpfs/data/fs71468/GenAI_para_runs/results/exp", 
                                    runs[run, "experiment"], "-run", runs[run, "run"],
                                    "/losses_tracker_1.csv"))
        result_2 <- read.csv(paste0("/gpfs/data/fs71468/GenAI_para_runs/results/exp", 
                                    runs[run, "experiment"], "-run", runs[run, "run"],
                                    "/losses_tracker_2.csv"))
        # tries since last valid result
        id_last_val <- max(which(!is.na(result_1$loss)))
        runs[run, "it_sinc_last_val_1"] <- nrow(result_1) - id_last_val
        id_last_val <- max(which(!is.na(result_2$loss)))
        runs[run, "it_sinc_last_val_2"] <- nrow(result_2) - id_last_val
        
        # quality criteria
        result <- rbind(result_1, result_2)
        best <- result[which.max(result$loss), ]
        best[, 4:8] <- round(best[, 4:8], 2)
        runs[run, 5:11] <- best[, c(1, 3:8)]
        runs[run, 14] <- max(result$iteration)
        runs[run, "num_evals"] <- nrow(result)
        runs[run, "num_valid_evals"] <- nrow(result[!is.na(result$loss), ])
    }, silent = TRUE)
    
    try({if(runs[run, "curr_iter"] > 19) restart <- FALSE}, silent = TRUE)
    try({if(runs[run, "curr_iter"] > 20 & runs[run, "status"] == "R") {
      restart <- FALSE
      system(paste0("scancel  ", runs[run, "id"]))
    }
    }, silent = TRUE)
    if(restart){
      if(!is.na(runs$id[run])) system(paste0("scancel ", runs[run, "id"]))
      bash_file <- paste0("/gpfs/data/fs71468/GenAI_para_runs/runs/exp",
                          runs[run, "experiment"], "-run", runs[run, "run"],
                          "/01_experiment_start.sh")
      system(paste0("sbatch ", bash_file))
    }
  }
  runs <- runs[, c(1:5, 14, 12:13, 6:11, 15:16)]
  
  # status report
  cat("\n\n*************************************************************************\n\n")
  print(runs)
  
  # start new node after 2.8 days
  if(difftime(Sys.time(), start_script, units = "secs")[[1]] > 60*60*24*2.8){
    keep_running <- FALSE
    setwd("/gpfs/data/fs71468/GenAI_para_runs/code/06_run_experiments")
    system("sbatch 02_experiment_master.sh")
  }
  Sys.sleep(60*60)
}






