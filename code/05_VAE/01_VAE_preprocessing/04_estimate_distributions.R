# Estimate TF distributions for GenAI_para
# GenAI_para-mHM project


# wd depending on system
if(Sys.info()[["nodename"]] == "cfgrammar"){
  main_path <- "/mnt/Data/Dropbox/Projekte/GenAI_para_for_HLSMs/"
} else {
  if(Sys.info()["sysname"] == "Linux"){
    main_path <- "/gpfs/data/fs71468/GenAI_para_runs/"
  }  
}
setwd(paste0(main_path, "data"))
code_dir <- paste0(main_path, "code")

code_files <- list.files(paste0(code_dir, "/02_utility_functions/09_CFG_implementation"), full.names = TRUE)
for(file in code_files) source(file)

setwd("07_VAE_data/03_prepared_TFs/prepared_functions")

# Create VSC skripts and folder for distribution estimation
lines <- readLines(paste0(code_dir, "/05_VAE/01_VAE_preprocessing/TEMPLATE_distribution_estimation.R"))
if(!dir.exists("../scripts/distribution_estimation_scripts")) dir.create("../scripts/distribution_estimation_scripts")

for (vae in c(1:5)){
  setwd(paste0(main_path, "data", "/07_VAE_data/03_prepared_TFs/prepared_functions"))
  folders <- list.files(pattern=paste0("vae_", vae))
  batch <- 0
  for(folder in folders){
    files <- list.files(folder, pattern = "_simplif")
    for(file in files){
      batch <- batch + 1
      lines[25] <- paste0('functions_for_dist <- "', folder, '/', file, '"')
      writeLines(lines,
                 con = paste0(
                   "../scripts/distribution_estimation_scripts/vae_", vae, 
                   "_distribution_estimation_batch",batch, ".R"))
    }
  }
  
  # Write bash file
  setwd("../scripts/distribution_estimation_scripts")
  fileConn <- file(paste0("vae_", vae, "distribution_estimation.sh"))
  writeLines(
    c("#!/bin/sh",
      "#SBATCH -J tf_distribution",
      "#SBATCH -N 1",
      "#SBATCH --qos=zen3_0512",
      "#SBATCH --partition=zen3_0512",
      paste0("#SBATCH --array=1-", batch, "%", min(c(40, batch))),
      "#SBATCH --mail-user=moritz.feigl@boku.ac.at",
      "#SBATCH --mail-type=BEGIN,END,FAIL",
      "source ~/env/spackenv",
      paste0("Rscript vae_", vae, "_distribution_estimation_batch${SLURM_ARRAY_TASK_ID}.R")),
    fileConn)
  close(fileConn)
  
  system(paste0("sbatch vae_", vae, "distribution_estimation.sh"))
}
