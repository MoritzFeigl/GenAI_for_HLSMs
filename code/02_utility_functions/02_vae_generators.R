# VAE generators for GenAI_para_for_HLSMs

# Generator functions --------------------------------------------------------------------
index_reconstructor <- function(pred_matrix){
  index <- integer(nrow(pred_matrix))
  for(i in 1:ncol(pred_matrix)) index[i] <- which.max(pred_matrix[, i])
  return(index)
}
index_sampler <- function(pred_matrix){
  index <- integer(nrow(pred_matrix))
  for(i in 1:ncol(pred_matrix)) index[i] <- sample(nrow(pred_matrix), 1, prob = pred_matrix[, i])
  return(index)
}

# helper function that makes the function evalutaion
tf_evaluation <- function(predicted_tf, variables){
  # paste variables=1 and function evaluation string
  var_def <- paste0(paste0(variables, "=1.0"), collaspe="
                    ")
  var_def <- paste0(var_def, collapse="")
  full_eval_fun <- paste0(var_def, "\n", predicted_tf)
                
  tf_eval <- try({
    eval(parse(text = paste('f_test <- function() {' ,  full_eval_fun , '}', sep='')))
    f_test()
  }, silent = TRUE)
  return(tf_eval)
}

variables <- c("dem", "aspect", "slope", "bd", "sand", "clay", # standard variables
               "lai", "mat_range", "map", "mat", # new variables
               "ThetaS", "KSat", "vGenu_n") # mHM variables

# main function for index prediction
generate_function_from_softmax <- function(vae, softmax_pred, data_dir, variables){
  
  # load vocabulary
  vocab <- read.csv(paste0(data_dir,"08_trained_vaes/vae", vae, "_vocabulary.csv"))
  vocab[vocab$X == 0, 2] <- ""
  vocab <- vocab[, 2]
  
  # transform log softmax to probabilities
  softmax_probabilities <- exp(softmax_pred)
  
  # get valid TF
  index_prediction <- index_reconstructor(softmax_probabilities)
  predicted_tf <- paste0(vocab[index_prediction], collapse = "")
  tf_eval <- tf_evaluation(predicted_tf, variables = variables)
  fail_count <- 0
  # random sample until a valid function is generated
  while(class(tf_eval) == "try-error" & fail_count < 2000){
    fail_count <- fail_count + 1
    index_prediction <- index_sampler(softmax_probabilities)
    predicted_tf <-  paste0(vocab[index_prediction], collapse = "")
    tf_eval <- tf_evaluation(predicted_tf, variables)
    if(is.null(tf_eval)) {
      tf_eval <- "try-error"
      class(tf_eval) <- "try-error"
    }
  }
  
  tf_eval <- tf_evaluation(predicted_tf, variables)
  if(class(tf_eval) == "try-error" | is.null(tf_eval)) return(NA)
  
  # fix ++ and -- if necessary
  for(i in 1:5) predicted_tf <- gsub("--", "-", predicted_tf, fixed = TRUE)
  for(i in 1:5) predicted_tf <- gsub("++", "-", predicted_tf, fixed = TRUE)
  
  return(predicted_tf)
}

# start vae in background
start_GenAI_para_vae <- function(run_dir, code_dir, vae_list, ncomplexes, 
                          vae_path = "/gpfs/data/fs71468/GenAI_para_runs/data/08_trained_vaes"){
  
  vae_pred_path <- paste0(run_dir, "/current_vae_predictions/")
  dir.create(vae_pred_path, showWarnings = FALSE)
  vaes <- unique(unlist(vae_list))
  
  # check if vae files for ncomplexes are available and create if not
  avail_decoder <- list.files(vae_path)
  for (vae in vaes){
    for (n in 1:ncomplexes){
      if(!(paste0("vae", vae, "_decoder_", n, ".pt") %in% avail_decoder)){
      file.copy(from = paste0(vae_path, "/vae", vae, "_decoder.pt"), 
                to = paste0(vae_path, "/vae", vae, "_decoder_", n, ".pt"))
      }
    }
  }
  for (setup in 1:ncomplexes){
    for (vae in vaes){
      
      # vae status files for python scripts
      status_file <- paste0(vae_pred_path, "vae", vae, "_setup", setup, "_start.txt")
      writeLines("end", status_file) # start key for vae script in background
      
      # vae generator python script
      pypath <- "export PYTHONPATH='/gpfs/data/fs71468/GenAI_para_for_HLSMs/GenAI_VAE'
           "
      script <- paste0("nohup /gpfs/data/fs71468/GenAI_para_for_HLSMs/miniconda3/envs/GenAI_VAE/bin/python3 ",
                       code_dir, "02_utility_functions/10_vae_generator.py ",
                       run_dir, " ")
      vae_cmd <- paste0(vae, " ", setup)
      
      # bash file
      fileConn <- file(paste0("vae", vae, "_", setup, ".sh"))
      writeLines(
        paste0(pypath, script, vae_cmd, " > vae", vae, "_", setup, ".out 2>&1 & 
                  exit", collapse = ""),
        fileConn)
      
      system(paste0(pypath, script, vae_cmd, " 1>/dev/null 2>/dev/null & 
                   exit", collapse = ""))
    }
  }
}



