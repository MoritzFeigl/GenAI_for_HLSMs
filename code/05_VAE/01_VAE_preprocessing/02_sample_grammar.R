# Sample Grammar for GenAI_para
# GenAI_para-ET project


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

# create Data and grammar directories
if(!dir.exists("07_VAE_data")) dir.create("07_VAE_data")
setwd("07_VAE_data")
if(!dir.exists("02_grammar_samples")) dir.create("02_grammar_samples")
setwd("02_grammar_samples")

# source utils
code_files <- list.files(paste0(code_dir, "/02_utility_functions/09_CFG_implementation"), full.names = TRUE)
for(file in code_files) source(file)

# Create grammar sampling scripts and corresponding bash files
lines <- readLines(paste0(code_dir, "/05_VAE/01_VAE_preprocessing/TEMPLATE_sample_grammar.R"))
n_batches <- 3
n_grammar_samples <- 40000000
seed_start <- sample(1000:10000, 1)
if(!dir.exists("scripts")) dir.create("scripts")
for(batch in 1:n_batches){
  lines[4] <- paste0('main_path <- "', main_path, '"')
  lines[24] <- paste0("funs <- grammar_sampler(n = ", 
                      ceiling(n_grammar_samples/n_batches), 
                      ", grammar = grammar, max_depth = 8,")
  lines[25] <- paste0("no_cores = 250, seed = ", seed_start+batch, ", save = TRUE,")
  lines[26] <- paste0('         file_name = paste0(main_path, "data/07_VAE_data/02_grammar_samples", sampled_grammar_batch_', batch,  '")')
  writeLines(lines,
             con = paste0(
               "scripts/sample_grammar_batch", batch, ".R"
             )
  )
}

# Write bash file
setwd("scripts")
fileConn <- file("sample_grammar_batches.sh")
writeLines(
  c("#!/bin/sh",
    "#SBATCH -J Grammar_sampling",
    "#SBATCH -N 1",
    "#SBATCH --qos=zen3_0512",
    "#SBATCH --partition=zen3_0512",
    paste0("#SBATCH --array=1-", n_batches, "%", min(c(10, n_batches))),
    "#SBATCH --mail-user=moritz.feigl@boku.ac.at",
    "#SBATCH --mail-type=BEGIN,END,FAIL",
    "source ~/env/spackenv",
    "Rscript sample_grammar_batch${SLURM_ARRAY_TASK_ID}.R"),
  fileConn)
close(fileConn)

# Run grammar sampling batches
system("sbatch sample_grammar_batches.sh")










