# generate developement data set
# GenAI_para-ET-mHM project

main_path <- "/gpfs/data/fs71468/GenAI_para_runs/"
code_dir <- paste0(main_path, "code")
setwd(paste0(main_path, "data/07_VAE_data/03_prepared_TFs/prepared_functions"))
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

no_cores <- 120

# define relevant dict of tokens
dict_from_list <- function(...){
  dictionary <- c(...)
  return(dictionary)
}
f <- c('sqrt', 'exp', 'log10', 'log', 
       'sin', 'sinh', 
       'cos', 'cosh', 
       'tan', 'tanh', 'atan',
       'abs')

base <- c("dem", "aspect_sin", "aspect_cos", "slope", "bd", "sand", "clay")
base_plus <- c("lai", "map", "mat", "mat_range")
mhm <- c("ThetaS", "KSat", "vGenu_n")
lai <- "lai"

vae_1 <- c(base)
vae_2 <- c(base, base_plus)
vae_3 <- c(base, mhm)
vae_4 <- c(base, base_plus, mhm)
vae_5 <- c(lai)
vae_version <- 1
operators <- c('+', '-', '*', '/', "^", "(", ")", ".")
numbers <- 0:9


dir.create(paste0(main_path, "data/07_VAE_data/04_vae_data"), showWarnings = FALSE)


for (vae in 1:5){
  vae_data_path <- paste0(main_path, "data/07_VAE_data/04_vae_data/vae_", vae)
  dir.create(vae_data_path, showWarnings = FALSE)
  count <- 0
  
  # vae specific dict of tokens
  var <- get(paste0("vae_", vae))
  dictionary <- dict_from_list(var, f, operators, numbers)
  put_in_spaces <- function(text){
    for(vocab in dictionary){
      text <- gsub(vocab, paste0(vocab, " "), text, fixed = TRUE)
      text <- gsub("sin h", paste0(vocab, " "), text, fixed = TRUE)
      text <- gsub("cos h", paste0(vocab, " "), text, fixed = TRUE)
      text <- gsub("tan h", paste0(vocab, " "), text, fixed = TRUE)
      text <- gsub("log 10", paste0(vocab, " "), text, fixed = TRUE)
      text <- gsub("mat _range", "mat_range", text, fixed = TRUE)
    }
    return(text)
  }
  
  
  batch_folders <- list.files(paste0(main_path, 
                                     "/data/07_VAE_data/03_prepared_TFs/prepared_functions"),
                              pattern = paste0("vae_", vae), full.names = TRUE)
  
  for (path in batch_folders){
    files <- list.files(path, pattern = "distribution", full.names = TRUE)
    
    for (file in files){
      count <- count + 1
      data <- fst::read_fst(file)
      data <- data[!(data$'function' %in% c("FALSE", "")), ]
      cat("file", count, " - nrow:", nrow(data), "\n")
      # put in spaces between functions
      prep_funs <- parallel::mclapply(X = data$'function',
                                      FUN = put_in_spaces,
                                      mc.cores = no_cores )
      data$'function'  <- unlist(prep_funs)
      
      # removing some unwanted TF types
      relicts <- c("Inf-", "Inf", "Inf+", "e-", "e+")
      for (remove in relicts){
        data <- data[!grepl(remove, data$'function', fixed = TRUE), ]
        
      }
      
      # remove TFs with no variance
      data <- data[(abs(data[, 2] - data[, 10])) > 0.1, ]
      # adapt functions that map to R- by multiplying with -1 and reverse quantiles
      neg_quantiles <- apply(data[, -1] < 0, 1, sum)
      data[which(neg_quantiles == 9), 1] <- paste0("- ( ", data[which(neg_quantiles == 9), 1], " )")
      data[which(neg_quantiles == 9), 2:10] <- data[which(neg_quantiles == 9), 10:2] * -1
      data <- data[neg_quantiles == 9 | neg_quantiles == 0, ]
      # remove functions that map out of [0, 1000]
      data <- data[data[, 10] <= 1000, ]
      
      file_name <- paste0(vae_data_path, "/vae_data_batch_", count, ".csv")
      write.csv(data, file_name, row.names = FALSE)
      rm(data)
    }
  }
}
