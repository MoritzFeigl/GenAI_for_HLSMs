# GenAI_para-ET-mHM project
# Create mHM compatible basin folder structures

library(ncdf4)
library(raster)
# Training/validation split for basins
if(Sys.info()[["nodename"]] == "cfgrammar"){
  setwd("/mnt/Data/Dropbox/Projekte/GenAI_para_for_HLSMs")
} else {
  if(Sys.info()["sysname"] == "Linux"){
    setwd("/gpfs/data/fs71468/GenAI_para_runs")
  } else {
    setwd("C:/Users/morit/Dropbox/Projekte/GenAI_para_for_HLSMs")
  }  
}

# Load basin infos
basins <- read.csv("data/06_training_validation_data/basin_split.csv")

# Training/validation folders
dir.create("data/06_training_validation_data/Training", showWarnings = FALSE)
dir.create("data/06_training_validation_data/Validation", showWarnings = FALSE)


for(basin in basins$Stat_ID){
  count <- which(basins$Stat_ID == basin)
  cat("Prepare basin", basin, paste0(count, "/", nrow(basins)), "\n")
  
  info <- basins[basins$Stat_ID == basin, 
                 c("Stat_ID", "Start_Date", "End_Date", "split")]
  if(info$split == "invalid") next
  for(i in 2:3) info[, i] <- as.Date(info[, i])
  
  # create folder in either Training or Validation
  basin_folder <- paste0("data/06_training_validation_data/", info$split, "/", info$Stat_ID)
  dir.create(basin_folder, showWarnings = FALSE)
  # config 
  # - mhm.nml, mhm_outputs.nml, mhm_parameter.nml, mpr.nml, mrm_outputs.nml, 
  #   run_single_basins.sh
  dir.create(paste0(basin_folder, "/config/"), showWarnings = FALSE)
  file.copy("data/05_mhm/master_mhm_config/mpr.nml",
            paste0(basin_folder, "/config/mpr.nml"), overwrite = TRUE)
  file.copy("data/05_mhm/master_mhm_config/mhm_outputs.nml",
            paste0(basin_folder, "/config/mhm_outputs.nml"), overwrite = TRUE)
  file.copy("data/05_mhm/master_mhm_config/mhm_parameter.nml",
            paste0(basin_folder, "/config/mhm_parameter.nml"), overwrite = TRUE)
  file.copy("data/05_mhm/master_mhm_config/mrm_outputs.nml",
            paste0(basin_folder, "/config/mrm_outputs.nml"), overwrite = TRUE)
  mhm_nml <- readLines("data/05_mhm/master_mhm_config/mhm.nml")
  mhm_nml <- gsub("basin_id", info$Stat_ID, mhm_nml)
  
  # set correct dates
  if(info$split == "Training"){
    
    # if end_date > after 2019, change it
    if(info$End_Date > "2020-01-01") info$End_Date <- as.Date("2019-12-31")
    
    mhm_nml[grep("dend = ", mhm_nml)] <- paste0("    eval_per(1)%dend = ", 
                                                format(info$End_Date, "%d"))
    mhm_nml[grep("dstart = ", mhm_nml)] <- paste0("    eval_per(1)%dstart = 01")
    mhm_nml[grep("mend = ", mhm_nml)] <- paste0("    eval_per(1)%mend = ", 
                                                format(info$End_Date, "%m"))
    mhm_nml[grep("mstart = ", mhm_nml)] <- paste0("    eval_per(1)%mstart = 01")
    mhm_nml[grep("yend = ", mhm_nml)] <- paste0("    eval_per(1)%yend = ", 
                                                format(info$End_Date, "%Y"))
    mhm_nml[grep("ystart = ", mhm_nml)] <- paste0("    eval_per(1)%ystart = 2014")
  }
  
  if(info$split == "Validation"){
    
    if(info$End_Date > "2013-12-31") info$End_Date <- as.Date("2013-12-31")
    # start either at 01.01.2000 or at the beginning of the time series
    if(as.POSIXct(info$Start_Date) > as.POSIXct("2005-01-01")){
      mhm_nml[grep("dstart = ", mhm_nml)] <- paste0("    eval_per(1)%dstart = ", 
                                                    format(info$Start_Date, "%d"))
      mhm_nml[grep("mstart = ", mhm_nml)] <- paste0("    eval_per(1)%mstart = ", 
                                                    format(info$Start_Date, "%m"))
      mhm_nml[grep("ystart = ", mhm_nml)] <- paste0("    eval_per(1)%ystart = ", 
                                                    format(info$Start_Date, "%Y"))
    } else {
      mhm_nml[grep("dstart = ", mhm_nml)] <- paste0("    eval_per(1)%dstart = 01")
      mhm_nml[grep("mstart = ", mhm_nml)] <- paste0("    eval_per(1)%mstart = 01")
      mhm_nml[grep("ystart = ", mhm_nml)] <- paste0("    eval_per(1)%ystart = 2005")
    }
    # end at 2013
    mhm_nml[grep("dend = ", mhm_nml)] <- paste0("    eval_per(1)%dend = ", 
                                                format(info$End_Date, "%d"))
    mhm_nml[grep("mend = ", mhm_nml)] <- paste0("    eval_per(1)%mend = ", 
                                                format(info$End_Date, "%m"))
    mhm_nml[grep("yend = ", mhm_nml)] <- paste0("    eval_per(1)%yend = ", 
                                                format(info$End_Date, "%Y"))
  }
  
  writeLines(mhm_nml, paste0(basin_folder, "/config/mhm.nml"))
  
  # forcings
  #   header.txt, pet.nc, pre.nc, tavg.nc
  dir.create(paste0(basin_folder, "/forcings/"), showWarnings = FALSE)
  for (forcing in c("pet.nc", "pre.nc", "tavg.nc", "header.txt")){
    file.copy(paste0("data/02_basin_data/sub_", info$Stat_ID, "/", forcing),
              paste0(basin_folder, "/forcings/", forcing), overwrite = TRUE)
  }
  
  # output --> empty
  dir.create(paste0(basin_folder, "/output"), showWarnings = FALSE)
  
  
  # static
  dir.create(paste0(basin_folder, "/static"), showWarnings = FALSE)
  
  
  # static/mpr
  #   spatial predictors
  dir.create(paste0(basin_folder, "/static/mpr"), showWarnings = FALSE)
  mpr_vars <- list.files(paste0("data/02_basin_data/sub_", info$Stat_ID, "/static/mpr/"))
  for (mpr_var in mpr_vars){
    file.copy(paste0("data/02_basin_data/sub_", info$Stat_ID, "/static/mpr/", mpr_var),
              paste0(basin_folder, "/static/mpr/", mpr_var), overwrite = TRUE)
  }
  ncin <- nc_open(paste0(basin_folder, "/static/mpr/land_cover.nc"), write = TRUE)
  lc_bounds <- ncvar_get(ncin,"land_cover_period_bnds")
  lc_bounds[2, 3] <- 2020
  ncvar_put(ncin,"land_cover_period_bnds", lc_bounds)
  nc_sync(ncin)
  nc_close(ncin)
  
  # LAI
  file.copy(paste0("data/02_basin_data/LAI/", info$split, "/lai_sub_", info$Stat_ID, ".nc"),
            paste0(basin_folder, "/static/mpr/lai.nc"), overwrite = TRUE)
  
  # LAI aggregated
  file.copy(paste0("data/02_basin_data/LAI_aggregated/", info$split, "/lai_sub_", info$Stat_ID, ".nc"),
            paste0(basin_folder, "/static/mpr/lai_agg.nc"), overwrite = TRUE)
  
  # MAT Range
  file.copy(paste0("data/02_basin_data/MAP_MAT_RANGE/", info$split, "/map_sub_", info$Stat_ID, ".nc"),
            paste0(basin_folder, "/static/mpr/map.nc"), overwrite = TRUE)
  file.copy(paste0("data/02_basin_data/MAP_MAT_RANGE/", info$split, "/mat_range_sub_", info$Stat_ID, ".nc"),
            paste0(basin_folder, "/static/mpr/mat_range.nc"), overwrite = TRUE) 
  file.copy(paste0("data/02_basin_data/MAP_MAT_RANGE/", info$split, "/mat_sub_", info$Stat_ID, ".nc"),
            paste0(basin_folder, "/static/mpr/mat.nc"), overwrite = TRUE)
  
  
  # static/routing: 
  #   [basin_id].txt (discharge file), dem.nc, facc.nc, fdir.nc, idgauges.nc
  dir.create(paste0(basin_folder, "/static/routing"), showWarnings = FALSE)
  routing_vars <- list.files(paste0("data/02_basin_data/sub_", info$Stat_ID, "/static/routing/"))
  for (routing_var in routing_vars){
    file.copy(paste0("data/02_basin_data/sub_", info$Stat_ID, "/static/routing/", routing_var),
              paste0(basin_folder, "/static/routing/", routing_var), overwrite = TRUE)
  }
  # check if new discharge is available and replace if it is the case
  discharge_file <- paste0(info$Stat_ID, ".txt")
  if(discharge_file %in% list.files("data/03_discharge/")){
    file.copy(paste0("data/03_discharge/", discharge_file),
              paste0(basin_folder, "/static/routing/", discharge_file), overwrite = TRUE)
  }
  # change discharge file format
  q_file <- readLines(paste0(basin_folder, "/static/routing/", discharge_file))
  second_line_id <- grep("nodata   -999.", q_file)
  q_file <- q_file[c(1, second_line_id:length(q_file))]
  writeLines(q_file, paste0(basin_folder, "/static/routing/", discharge_file))
  
  # ET
  dir.create(paste0(basin_folder, "/ET"), showWarnings = FALSE)
  file.copy(paste0("data/02_basin_data/ET/", info$split, "/et_sub_", basin, ".nc"),
            paste0(basin_folder, "/ET/et.nc"), overwrite = TRUE)
  
}

# Fix lai names
setwd("data/06_training_validation_data")
train_basins <- list.files("Training", full.names = TRUE)
val_basins <- list.files("Validation", full.names = TRUE)
basins <- c(train_basins, val_basins)

for(basin in basins){
  path <- paste0(basin, "/static/mpr/lai.nc")
  ncin <- nc_open(path, write=TRUE)
  try({ncin <- ncvar_rename(ncin, "lai_class", "lai", verbose=FALSE )})
  nc_sync(ncin)
  nc_close(ncin)
  
  path <- paste0(basin, "/static/mpr/lai_agg.nc")
  ncin <- nc_open(path, write=TRUE)
  try({ncin <- ncvar_rename(ncin, "lai_class", "lai_agg", verbose=FALSE )})
  nc_sync(ncin)
  nc_close(ncin)
}



# fix missing values in 9315014 and validation basins
val_basins <- paste0("Validation/", list.files("Validation"))
fix_basins <- c("Training/9315014", val_basins)

for(basin in fix_basins){
  for (var in c("mat", "mat_range", "map")){
    
    if(!file.exists(paste0(basin, "/static/mpr/", var, " - kopie.nc"))){
      # copy var
      file.copy(paste0(basin, "/static/mpr/", var, ".nc"), 
                paste0(basin, "/static/mpr/", var, " - kopie.nc"))
      
      # copy slope
      path <- paste0(basin, "/static/mpr/", var, ".nc")
      file.copy(paste0(basin, "/static/mpr/slope.nc"), 
                path, overwrite = TRUE)
      
      # rename slope to var
      ncin <- nc_open(path, write=TRUE)
      try({ncin <- ncvar_rename(ncin, "slope", var, verbose=FALSE )})
      nc_sync(ncin)
      nc_close(ncin)
      
      # load original var values
      or_var <- nc_open(paste0(basin, "/static/mpr/", var, " - kopie.nc"), write=TRUE)
      old_var <- ncvar_get(or_var, var)
      
      # put in var values in renamed slope nc
      ncin <- nc_open(path, write=TRUE)
      new_var <- ncvar_get(ncin, var)
      relevant_ids <- which(!is.na(new_var), arr.ind = TRUE)
      new_var[relevant_ids] <- old_var[relevant_ids]
      w <- matrix(1, 11, 11)
      a <- raster(old_var)
      x <- focal(a, w, mean, na.rm=TRUE, NAonly=TRUE, pad=TRUE)
      xx <- focal(x, w, mean, na.rm=TRUE, NAonly=TRUE, pad=TRUE)
      new_var[relevant_ids] <- xx[relevant_ids]
      ncvar_put(ncin, var, new_var)
      nc_sync(ncin)
      nc_close(ncin)
    }
  }
}



