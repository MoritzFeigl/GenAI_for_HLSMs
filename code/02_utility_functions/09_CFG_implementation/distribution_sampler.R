#' distribution_sampler
#'
#' @param functions
#' @param variable_df
#' @param scaling_bounds
#' @param file_name
#' @param no_cores
#'
#' @return
#' @export
#'
#' @examples
distribution_sampler <- function(functions, variable_df,
                                 scaling_bounds = NULL,
                                 file_name = NULL,
                                 no_cores = NULL,
                                 batches = 20){
  
  options(scipen = 999)
  
  # check for NAs
  if(sum(is.na(variable_df)) != 0) {
    stop(paste("NA values detected in variable_df!",
               "\nPlease provide a data frame with either imputated NA values",
               "or remove rows with NAs."))
  }
  
  if(scaling_bounds != FALSE){
    if(sum(names(variable_df) %in% names(scaling_bounds)) != ncol(variable_df)){
      stop(paste("Not all variables of variable_df are included in the provided scaling_bounds",
                 "\nPlease provide scaling bounds that fit variable_df, or check column names."))
    }
    # scale variables
    add_bounds <- names(variable_df)[!(names(variable_df) %in% names(scaling_bounds))]
    for(variable in add_bounds){
      scaling_bounds[[variable]] <- c(min(variable_df[, variable]),
                                      max(variable_df[, variable]))
    }
    for(variable in names(variable_df)){
      variable_df[, variable] <- .range01(variable_df[, variable],
                                          scaling_bounds[[variable]])
      variable_df[, variable][variable_df[, variable] < 0.01] <- 0.01
    }
  }
  if(is.data.frame(functions)){
    if(ncol(functions) == 1){
      functions <- as.character(functions[, 1])
    } else stop("functions must be either a character vector or a 1 column data frame with function strings.")
  } else {
    if(!is.character(functions)){
      stop("functions must be either a character vector or a 1 column data frame with function strings.")
    }
  }
  
  # define file_name if not specified
  if(is.null(file_name)){
    file_name <- paste0("variable_input_grammar", "-",
                        format(Sys.time(), "%d-%m-%Y-%H%M"))
  }
  # create folder
  file_name <- gsub(".fst", "", file_name)
  cat("Results will be saved in", paste0(file_name, "_batches [0-", batches, "]"))
  #if(!dir.exists(file_name))  dir.create(file_name)
  
  
  # split in batches
  if(batches > 1){
    batch_cuts <- cut(1:length(functions), breaks = batches)
    functions_batches <- split(functions, batch_cuts)
  } else {
    functions_batches <- list(functions)
  }
  
  for(batch in 1:batches){
    
    # parallel computation
    parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
    if(Sys.info()[["sysname"]] == "Windows"){
      cl <- parallel::makeCluster(no_cores)
      parallel::clusterExport(cl, list(".evaluate_function_from_string",
                                       ".distribution_values_from_fun",
                                       "variable_df"
      ),
      envir = environment())
      output_list <- pbapply::pblapply(functions_batches[[batch]],
                                       FUN = .distribution_values_from_fun,
                                       variable_df = variable_df,
                                       cl = cl)
      parallel::stopCluster(cl)
      output <- do.call(rbind, output_list)
      rownames(output) <- NULL
    } else {
      output_list <- parallel::mclapply(functions_batches[[batch]],
                                        FUN = .distribution_values_from_fun,
                                        variable_df = variable_df,
                                        mc.cores = no_cores )
      output <- do.call(rbind, output_list)
      rownames(output) <- NULL
    }
    
    # Remove functions that are NA/Inf
    for (col in names(output)){
      output <- output[!is.na(output[, col]), ]
      output <- output[!is.infinite(output[, col]), ]
    }
    
    # save batch
    batch_name <- paste0(file_name, "_batch", batch)
    fst::write_fst(x = output, path = batch_name, compress = 100)
    rm(output)
  }
}
