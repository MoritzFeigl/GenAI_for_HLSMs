# GenAI_para-ET-mHM project
# analyze validation


# load geos, proj and gdal for raster
system("spack load /54rh3fo
spack load /da7sz6h
spack load /4b4eze3
export LD_LIBRARY_PATH=$LIBRARY_PATH")
# load libraries
library(raster)
library(ggplot2)
library(rasterVis)
library(magrittr)
library(RColorBrewer)
library(viridis)  # better colors for everyone
#library(ggthemes) # theme_map()
library(ncdf4)
library(ggplot2)
library(gridExtra)
#library(stars)
# 1. paths & colors ----------------------------------------------------------------------
for(split in c("Validation", "Training", "Validation_2km", "Validation_1km")){
  
  run_name <- split
  main_path <- "/gpfs/data/fs71468/GenAI_para_runs/"
  run_dir <- paste0(main_path, "runs/", run_name)
  result_dir <- paste0(main_path, "results/", run_name)
  code_dir <- paste0(main_path, "code/")
  data_dir <- paste0(main_path, "data/")
  setup <- n_parallel_setups <- 1
  dir.create(run_dir, showWarnings = FALSE)
  dir.create(result_dir, showWarnings = FALSE)
  source(paste0(code_dir, "02_utility_functions/01_GenAI_para_utility_functions.R"))
  training_results_file <- paste0(main_path, "results/training_results_2024-01-16.csv")
  
  # read netcdfs ---------------------------------------------------------------------------
  training_results <- read.csv(training_results_file)
  mhm_results <- unique(training_results[training_results$experiment == 0, 
                                         c("run", "mean_loss", "time")])
  mhm_run <- mhm_results[mhm_results$mean_loss == max(mhm_results$mean_loss), "run"]
  
  
  for (experiment in c(1, 2)){
    # get best experiment/run combinations
    experiment_name <- c("1" = "GenAI_para",
                         "2" = "GenAI_para + SP")[as.character(experiment)]
    
    exp_results <- unique(training_results[training_results$experiment == experiment, 
                                           c("run", "mean_loss", "time")])
    run <- exp_results[exp_results$mean_loss == max(exp_results$mean_loss), "run"]
    
    # mHM vs GenAI_para parameter fields
    exp_path <- paste0(result_dir, "/exp", experiment, "-run", run)
    mhm_path <- paste0(result_dir, paste0("/exp0-run", mhm_run))
    para_files <- list.files(path = exp_path, 
                             pattern = "parameters.nc")
    basins <- sapply(para_files, function(x) gsub("_parameters.nc", "", x))
    plot_path <- paste0(main_path, "/analysis_results/", tolower(split), 
                        "_results/parameter_plots/") 
    dir.create(plot_path, recursive = TRUE)
    
    # function to cretae df from raster
    raster_to_df <- function(raster, flip_vertical = FALSE, flip_horizontal = FALSE) {
      if (flip_vertical) {
        raster <- flip(raster, direction = 'y')  # Flip vertically
      }
      if (flip_horizontal) {
        raster <- flip(raster, direction = 'x')  # Flip horizontally
      }
      df <- as.data.frame(raster, xy = TRUE)
      names(df) <- c("lon", "lat", "value")
      return(df)
    }
    
    for(basin in basins){
      
        basin_path <- paste0(plot_path, basin)
        dir.create(basin_path, recursive = TRUE)
        
        
        # Parameter plots ----------------------------------------------------------------------
        parameters <- c("KSat_till", 
                        "FieldCap_till",
                        "L1_FieldCap", 
                        "fRoots_temp",
                        "L1_fRoots",
                        "ThetaS_till", 
                        "L1_SatSoilMoisture", 
                        "L1_Max_Canopy_Intercept")
        
        for(parameter in parameters){
          try({
          # GenAI_para raster
          basin_nc <- nc_open(paste0(exp_path, "/", basin, "_parameters.nc"))
          
          # Extract parameter data
          param_data <-  ncvar_get(basin_nc, parameter)
          if(parameter == "L1_Max_Canopy_Intercept"){
            # in case of max canopa get july data and remove values outside the basin
            GenAI_para_raster_nc <- raster(t(param_data[, , 7]))
            
            satsoilm_data <-  ncvar_get(basin_nc, "L1_SatSoilMoisture")
            satsoilm_nc <- raster(t(satsoilm_data[1, , , 1]))
            values(GenAI_para_raster_nc)[is.na(values(satsoilm_nc))] <- NA
          } else {
            GenAI_para_raster_nc <- raster(t(param_data[1, , , 1]))
          }
          
          #########
          GenAI_para_raster_nc <- raster(param_data[1, , , 1])
          
          
          
          
          lat <-  ncvar_get(basin_nc, "lat_bnds")
          lon <-  ncvar_get(basin_nc, "lon_bnds")
          resx <- 100
          resy <- 100
          xmn <- min(lon)# + 0.5 * resx
          xmx <- max(lon)# - 0.5 * resx
          ymn <- min(lat)# + 0.5 * resy
          ymx <- max(lat)# - 0.5 * resy
          

          # Create a template raster (done once outside the loop)
          r1 <- terra::rast(terra::ext(xmn, xmx, ymn, ymx), resolution=c(resx, resy))
          crs(r1) <- "EPSG:31467"
          
          # Fill temporary raster with TWSA values for the current time step
          r1[] <- as.numeric(t(param_data[1, , , 1]))
          
          
          
          # Save the plots side by side
          png(paste0(basin_path, "/", parameter, "_exp", experiment, "-run", run, "_TEST.png"),
              width = 25, height = 10, units = "in", res = 150)
          plot(r1)
          dev.off()
          
          
          ##########
          lat_data <-  ncvar_get(basin_nc, "lat_bnds")
          lon_data <-  ncvar_get(basin_nc, "lon_bnds")
          extent(GenAI_para_raster_nc) <- c(min(lat_data), max(lat_data), min(lon_data), max(lon_data))
          crs(GenAI_para_raster_nc) <- CRS("+init=epsg:31467")
          nc_close(basin_nc)
          
          # mHM raster
          basin_nc <- nc_open(paste0(mhm_path, "/", basin, "_parameters.nc"))
          param_data <-  ncvar_get(basin_nc, parameter)
          if(parameter == "L1_Max_Canopy_Intercept"){
            mhm_raster_nc <- raster(t(param_data[, , 7]))
            values(mhm_raster_nc)[is.na(values(satsoilm_nc))] <- NA
            
          } else {
            mhm_raster_nc <- raster(t(param_data[1, , , 1]))
          }
          lat_data <-  ncvar_get(basin_nc, "lat_bnds")
          lon_data <-  ncvar_get(basin_nc, "lon_bnds")
          extent(mhm_raster_nc) <- c(min(lat_data), max(lat_data), min(lon_data), max(lon_data))
          nc_close(basin_nc)
          
          # save Raster
          GenAI_para_raster_flipped <- flip(GenAI_para_raster_nc, direction = 'y')
          writeRaster(GenAI_para_raster_flipped, 
                      filename = paste0(basin_path, "/", parameter, "_exp", experiment, 
                                        "-run", run, ".tif"), 
                      format = "GTiff", overwrite=TRUE)
          
          # Convert rasters to data frames, adjust flip as necessary
          GenAI_para_df <- raster_to_df(GenAI_para_raster_nc, flip_vertical = TRUE)
          
          mhm_df <- raster_to_df(mhm_raster_nc, flip_vertical = TRUE)
          
          # Create ggplot objects
          if(parameter == "KSat_till" & experiment == 2){
            
          p1 <- ggplot(GenAI_para_df, aes(x = lon, y = lat, fill = value)) +
            geom_raster(na.rm = TRUE) +
            scale_fill_viridis(name = parameter, na.value = NA, option = "D", trans = "log") +
            theme_minimal() +
            labs(title = experiment_name)
          } else {
              p1 <- ggplot(GenAI_para_df, aes(x = lon, y = lat, fill = value)) +
                geom_raster(na.rm = TRUE) +
                scale_fill_viridis(name = parameter, na.value = NA, option = "D") +
                theme_minimal() +
                labs(title = experiment_name)
            }
          
          p2 <- ggplot(mhm_df, aes(x = lon, y = lat, fill = value)) +
            geom_raster(na.rm = TRUE) +
            scale_fill_viridis(name = parameter, na.value = NA, option = "D") +  
            theme_minimal() +
            labs(title = "mHM standard functions")
          
          # Save the plots side by side
          png(paste0(basin_path, "/", parameter, "_exp", experiment, "-run", run, "_vs_mHM.png"),
              width = 25, height = 10, units = "in", res = 150)
          grid.arrange(p1, p2, ncol = 2)
          dev.off()
          })
        }
        
        # ET plot ----------------------------------------------------------------------------
        # get lat lon data
        try({
        basin_nc <- nc_open(paste0(exp_path, "/", basin, "_parameters.nc"))
        lat_data <-  ncvar_get(basin_nc, "lat_bnds")
        lon_data <-  ncvar_get(basin_nc, "lon_bnds")
        nc_close(basin_nc)
        
        # get all July ET patterns and aggregate
        get_mean_month_et <- function(data, month){
          year_months <- NULL
          years <- dim(data)[3]/12
          for(i in 0:(years-1)) year_months[i+1] <- i*12 +7
          return(apply(param_data[, , year_months], c(1, 2), mean))
        }
        
        month <- 7
        
        # GenAI_para raster
        basin_nc <- nc_open(paste0(exp_path, "/", basin, "_Fluxes_States.nc"))
        param_data <-  ncvar_get(basin_nc, "aET")
        GenAI_para_raster_nc <- raster(t(get_mean_month_et(param_data, month)))
        extent(GenAI_para_raster_nc) <- c(min(lat_data), max(lat_data), min(lon_data), max(lon_data))
        crs(GenAI_para_raster_nc) <- CRS("+init=epsg:31467")
        nc_close(basin_nc)
        
        # mHM raster
        basin_nc <- nc_open(paste0(mhm_path, "/", basin, "_Fluxes_States.nc"))
        param_data <-  ncvar_get(basin_nc, "aET")
        mhm_raster_nc <- raster(t(get_mean_month_et(param_data, month)))
        extent(mhm_raster_nc) <- c(min(lat_data), max(lat_data), min(lon_data), max(lon_data))
        crs(mhm_raster_nc) <- CRS("+init=epsg:31467")
        nc_close(basin_nc)
        
        # ALEXIS raster
        alexi_split <- strsplit(split, "_")[[1]][1]
        basin_nc <- nc_open(paste0(main_path, "data/02_basin_data/ET/", alexi_split, "/et_sub_", 
                                   basin, ".nc"))
        param_data <-  ncvar_get(basin_nc, "__xarray_dataarray_variable__")
        ALEXIS_raster_nc <- raster(t(get_mean_month_et(param_data, month)))
        values(ALEXIS_raster_nc)[is.na(values(GenAI_para_raster_nc))] <- NA
        extent(ALEXIS_raster_nc) <- c(min(lat_data), max(lat_data), min(lon_data), max(lon_data))
        crs(ALEXIS_raster_nc) <- CRS("+init=epsg:31467")
        nc_close(basin_nc)
        
        # Convert rasters to data frames, adjust flip as necessary
        GenAI_para_df <- raster_to_df(GenAI_para_raster_nc, flip_vertical = TRUE)
        mhm_df <- raster_to_df(mhm_raster_nc, flip_vertical = TRUE)
        ALEXIS_raster_df <- raster_to_df(ALEXIS_raster_nc, flip_vertical = TRUE)
        
        # Create ggplot objects
        p1 <- ggplot(GenAI_para_df, aes(x = lon, y = lat, fill = value)) +
          geom_raster(na.rm = TRUE) +
          scale_fill_viridis(name = "aET", na.value = NA, option = "D") +
          theme_minimal() +
          labs(title = experiment_name)
        
        p2 <- ggplot(mhm_df, aes(x = lon, y = lat, fill = value)) +
          geom_raster(na.rm = TRUE) +
          scale_fill_viridis(name = "aET", na.value = NA, option = "D") +  
          theme_minimal() +
          labs(title = "mHM standard functions")
        
        p3 <- ggplot(ALEXIS_raster_df, aes(x = lon, y = lat, fill = value)) +
          geom_raster(na.rm = TRUE) +
          scale_fill_viridis(name = "aET", na.value = NA, option = "D") +  
          theme_minimal() +
          labs(title = "ALEXI")
        
        
        # Save the plots side by side
        png(paste0(basin_path, "/aET_exp", experiment,"-run", run, "_vs_mHM.png"), 
            width = 30, height = 10, units = "in", res = 150)
        grid.arrange(p1, p2, p3, ncol = 3)
        dev.off()
      })
    }
  }
}

