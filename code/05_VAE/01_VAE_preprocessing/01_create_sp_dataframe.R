# GenAI_para-mHM-ET project
# Create spatial predictors data frame from netcdf


if("RNetCDF" %in% rownames(installed.packages()) == FALSE) {
  install.packages("RNetCDF")
}
if("feather" %in% rownames(installed.packages()) == FALSE) {
  install.packages("feather")
}
library(RNetCDF)

if(Sys.info()[["nodename"]] == "cfgrammar"){
  setwd("/mnt/Data/Dropbox/Projekte/GenAI_para_for_HLSMs/data")
} else {
  if(Sys.info()["sysname"] == "Linux"){
    setwd("/gpfs/data/fs71468/GenAI_para_runs/data")
  }  
  if(Sys.info()["sysname"] == "Windows"){
    setwd("C:/Users/morit/Dropbox/Projekte/GenAI_para_for_HLSMs/data")
  }  
}
if(!dir.exists("/07_VAE_data")) dir.create("07_VAE_data")
if(!dir.exists("/07_VAE_data/01_sp_data")) dir.create("07_VAE_data/01_sp_data")
basins <- list.files("06_training_validation_data/Training")
for (basin in basins){
  # load mpr input netcdf, aggregate layers and combine in data.frame
  cat("******", basin, "******\n")
  variables <- c("mpr/bd", "mpr/sand", "mpr/clay", "mpr/slope", "mpr/aspect", 
                 "routing/dem", "mpr/lai", "mpr/map", "mpr/mat", "mpr/mat_range")
  for(variable in variables){
    cat(variable, "\n")
    nc <- open.nc(
      paste0("06_training_validation_data/Training/",
             basin, "/static/", variable, ".nc")
    )
    if(variable == "mpr/lai"){
      var <- var.get.nc(nc, "lai_class")
    } else {
    var <- var.get.nc(nc, unlist(strsplit(variable, "/"))[2])
    }
    if(length(dim(var)) == 3) var <- apply(var, c(2, 3), function(x) mean(x, na.rm = TRUE))
    var_vector <- as.data.frame(as.vector(var))
    names(var_vector) <- unlist(strsplit(variable, "/"))[2]
    if(!exists("spatial_predictors")) {
      spatial_predictors <- var_vector
    } else {
      spatial_predictors <- cbind(spatial_predictors, var_vector)
    }
    close.nc(nc)
  }
  
  # mHM specific variables
  variables <- c("KSat_till", "KSat_notill", "vGenu_n_till", "vGenu_n_notill", "ThetaS_till", "ThetaS_notill")
  nc <- open.nc(paste0("../runs/03_standard_mhm_para_run/Training/", basin, "/output/mHM_parameters.nc"))
  for(variable in variables){
    cat(variable, "\n")
    var <- var.get.nc(nc, variable)
    if(length(dim(var)) > 2) var2 <- apply(var, c(2, 3), function(x) mean(x, na.rm = TRUE))
    var_vector <- as.data.frame(as.vector(var))
    names(var_vector) <- variable
    spatial_predictors <- cbind(spatial_predictors, var_vector)
  }
  nas <- which(is.na(spatial_predictors), arr.ind = TRUE)
  if(length(nas) != 0){ 
    spatial_predictors <- spatial_predictors[-nas[, 1], ]
  }
  feather::write_feather(spatial_predictors,
                         paste0("07_VAE_data/01_sp_data/",
                                basin, "_spatial_predictors.feather"))
  rm(spatial_predictors)
}

all_sp_files <- list.files("07_VAE_data/01_sp_data", full.names = TRUE, pattern = ".feather")
all_sp_files <- all_sp_files[all_sp_files != "spatial_predictors.feather"]
cat("Aggregating spatial predictor data frames\n")
for(basin in basins){
  cat(basin, "\n")
  sp_curr <- feather::read_feather(all_sp_files[grep(basin, all_sp_files)])
  sp_curr$basin <- basin
  if(!exists("spatial_predictors")) {
    spatial_predictors <- sp_curr
  } else {
    spatial_predictors <- rbind(spatial_predictors, sp_curr)
  }
  rm(sp_curr)
}

feather::write_feather(spatial_predictors,"07_VAE_data/01_sp_data/spatial_predictors.feather")


# Create sampled sp data.frame
spatial_predictors <- feather::read_feather("07_VAE_data/01_sp_data/spatial_predictors.feather")
for(basin in basins){
  basin_data <- spatial_predictors[spatial_predictors$basin == basin, ]
  # sample always 0.25 % of the basin cells to reach approx. 10 000 sampled data points
  samples <- ceiling(nrow(basin_data)*0.0025)
  if(exists("sp_df_sampled")){
    sp_df_sampled <- rbind(sp_df_sampled,
                           basin_data[sample(nrow(basin_data), samples), ])
  } else {
    sp_df_sampled <- basin_data[sample(nrow(basin_data), samples), ]
  }
}

feather::write_feather(sp_df_sampled,
                       "07_VAE_data/01_sp_data/spatial_predictors_sampled.feather")

# remove all individual basin files
for(basin in basins){
  unlink(paste0("07_VAE_data/01_sp_data/", basin, "_spatial_predictors.feather"))
}
