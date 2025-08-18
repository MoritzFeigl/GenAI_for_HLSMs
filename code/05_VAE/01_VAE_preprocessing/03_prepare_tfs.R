# Prepare transfer functions for GenAI_para
# GenAI_para-ET-mHM project


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

# get Grammar samples
setwd("07_VAE_data/02_grammar_samples")
n_batches <- length(list.files(pattern=".fst"))
for(batch in 1:n_batches){
  sampled_batch <- fst::read_fst(paste0("sampled_grammar_batch_", batch, ".fst"))
  if(exists("grammar_samples")){
    grammar_samples <- rbind(grammar_samples, sampled_batch)
  } else grammar_samples <- sampled_batch
}

grammar_samples <- data.frame(functions = unique(grammar_samples$functions))
fst::write_fst(grammar_samples, "sampled_grammar.fst")
file.remove(paste0("sampled_grammar_batch_", 1:n_batches, ".fst"))

setwd("..")
if(!dir.exists("03_prepared_TFs")) dir.create("03_prepared_TFs")
setwd("03_prepared_TFs")

# Split sampled grammar into multiple batches
create_fst_batch_files <- function(files, number_of_batches_per_file){
  # Script to split available fst files in smaller batches for distributed computing
  # split all files
  count <- 0
  file_names <- gsub(".fst", "",
                     tail(unlist(strsplit(files[1], split = "/")), 1), fixed = TRUE)
  for(k in seq_along(files)){
    if(!dir.exists(paste0(file_names[k], "_batches"))){
      dir.create(paste0(file_names[k], "_batches"))
    }
    functions_para <- fst::read_fst(files[k])
    batch_cuts <- cut(1:nrow(functions_para),
                      breaks = number_of_batches_per_file)
    functions_batches <- split.data.frame(functions_para, batch_cuts)
    for(i in 1:number_of_batches_per_file){
      count <- count + 1
      fst::write.fst(
        functions_batches[[i]],
        path = paste0(file_names[k], "_batches/", file_names[k], "_batch",
                      count, ".fst"))
    }
  }
}

create_fst_batch_files("../02_grammar_samples/sampled_grammar.fst", 2)


# Create grammar sampling scripts and corresponding bash files
lines <- readLines(paste0(code_dir, "/05_VAE/01_VAE_preprocessing/TEMPLATE_prepare_tfs.R"))
n_batches <- length(list.files("sampled_grammar_batches"))
seed_start <- sample(1000:10000, 1)
if(!dir.exists("scripts")) dir.create("scripts")


# VAE versions
# 1: base
# 2: base & base_plus
# 3: base & mhm
# 4: base & base_plus & mhm
# 5: lai
for(version in 1:5){
  for(batch in 1:n_batches){
    # KSat variable inputs and simplify
    lines[29] <- paste0("vae_version <- ", version)
    lines[36] <- paste0('variable_input(functions = "', "sampled_grammar_batches/",
                        "sampled_grammar_batch",
                        batch, '.fst",')
    lines[40] <- paste0("n_iter = 20, no_cores = 110, seed = ", seed_start+batch, ",")
    lines[41] <- paste0('         file_name = "prepared_functions/vae_', 
                        version, '_functions_batch_', batch,  '")')
    lines[42] <- paste0('setwd("prepared_functions/vae_', 
                        version, '_functions_batch_', batch,  '")')
    lines[43] <- paste0('simplify_functions(files = paste0("vae_', 
                        version, '_functions_batch_', batch, '_batch", 1:10, ".fst"),')
    writeLines(lines, con = paste0("scripts/vae_", version, "_batch_", batch, ".R"))
  }
}

# create folder for prepared TFs
if(!dir.exists("prepared_functions")) dir.create("prepared_functions")

# Write bash file
setwd("scripts")
for(version in 1:5){
fileConn <- file(paste0("prepare_tfs_vae_", version, ".sh"))
  writeLines(
    c("#!/bin/sh",
      paste0("#SBATCH -J vae_", version, "_tfs"),
      "#SBATCH -N 1",
      "#SBATCH --qos=zen3_0512",
      "#SBATCH --partition=zen3_0512",
      paste0("#SBATCH --array=1-", n_batches, "%", min(c(20, n_batches))),
      "#SBATCH --mail-user=moritz.feigl@boku.ac.at",
      "#SBATCH --mail-type=BEGIN,END,FAIL",
      "source ~/env/spackenv",
      paste0("Rscript vae_", version, "_batch_${SLURM_ARRAY_TASK_ID}.R")),
    fileConn)
  close(fileConn)
}

# Run grammar sampling batches
for(version in 1:5){
  system(paste0("sbatch prepare_tfs_vae_", version, ".sh"))
}



