
# Print all TFs for specific experiment
exp <- 3
for(run in 1:5){
  setwd(paste0("/gpfs/data/fs71468/GenAI_para_runs/results/exp", exp,"-run", run))
  tfs <- tail(read.csv("best_TFs_tracker_1.csv"), 1)
  
  cat("\n", "--------------", paste0("exp", exp,"-run", run), "--------------", "\n")
  for(i in 3:11) cat(names(tfs)[i], ":", as.character(tfs[i]), "\n")
}






# check number of valid experiments
exp <- 1
run <- 1
setwd(paste0("/gpfs/data/fs71468/GenAI_para_runs/results/exp", exp, "-run", run))

results <- NULL
for(setup in 1:2){
  loss <- read.csv(paste0("losses_tracker_", setup, ".csv"))
  tfs <- read.csv(paste0("TFs_tracker_", setup, ".csv"))
  results <- rbind(results, merge(loss, tfs, by="time"))
}

rm_row <- unique(which(is.na(results[, c("Ksat", "FieldCap", "ThetaS_1", "ThetaS_2", "fRoots_1", 
                        "fRoots_2", "L1_Max_Canopy_Intercept")]), arr.ind = TRUE)[, 1])
check <- results[-rm_row, ]
check <- check[is.na(check$loss), ]
check$time <- as.date(check$time)
sub_check <- check[check$time > "2023-08-07", ]


valid_results <- results[!is.na(results$loss), ]
valid_results[order(valid_results$loss, decreasing = TRUE), 4:8]


# check weird outcomes
one <- check[1, ]
one$time

parameters <- NULL
for(setup in 1:2){
  para <- read.csv(paste0("parameters_tracker_", setup, ".csv"))
  parameters <- rbind(parameters, para)
}
# check KSat prediction
one_para <- parameters[parameters$time == one$time, ]
point <- one_para[grep("GenAI_para", names(one_para))]


