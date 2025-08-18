setwd("C:/Users/morit/Dropbox/Projekte/GenAI_para_for_HLSMs/data/10_compare_data")

library(ncdf4)
library(raster)
library(rnaturalearth)
library(sf)

# Germany mask
ger_sf <- ne_countries(country = "Germany", returnclass = "sf")

# helper to get Date vector from a GLEAM file
get_dates <- function(ncfile){
  nc    <- nc_open(ncfile)
  tvec  <- ncvar_get(nc, "time")
  tunits <- ncatt_get(nc, "time", "units")$value
  origin <- sub("days since ", "", tunits)
  dates <- as.Date(tvec, origin = origin)
  nc_close(nc)
  dates
}

years     <- 2005:2013
max_julys <- vector("list", length(years))

for(i in seq_along(years)){
  print(i)
  yr     <- years[i]
  file   <- sprintf("Ei_%d_GLEAM_v4.2a.nc", yr)
  
  # 1) load as RasterBrick and attach real dates
  rb     <- brick(file, varname = "Ei")
  dates  <- get_dates(file)
  rb     <- setZ(rb, dates, name = "time")
  names(rb) <- as.character(dates)
  
  # 2) crop & mask to Germany
  rc     <- crop(rb, ger_sf)
  rm     <- mask(rc, ger_sf)
  
  # 3) select July layers
  juli   <- which(format(getZ(rm), "%m") == "07")
  
  # 4) compute the max over July days
  max_julys[[i]] <- calc(rm[[juli]], fun = max, na.rm = TRUE)
}

# stack all yearly maxima and take the mean
s      <- stack(max_julys)
mean_july_max <- calc(s, fun = max, na.rm = TRUE)


# 8) Export raster to GeoTIFF
crs(mean_july_max) <- "EPSG:4326"
output_path <- "C:/Users/morit/Dropbox/Projekte/GenAI_para_for_HLSMs/analysis_results/compare_tiffs/"
dir.create(output_path, showWarnings = FALSE)
output_file <- "max_intercept_Germany_GLEAM.tiff"
writeRaster(mean_july_max, filename=paste0(output_path, output_file), overwrite=TRUE, NAflag=-9999)


# also mean max
mean_july_mean_max <- calc(s, fun = mean, na.rm = TRUE)


# 8) Export raster to GeoTIFF
crs(mean_july_mean_max) <- "EPSG:31467"
output_path <- "C:/Users/morit/Dropbox/Projekte/GenAI_para_for_HLSMs/analysis_results/compare_tiffs/"
dir.create(output_path, showWarnings = FALSE)
output_file <- "mean_max_intercept_Germany_GLEAM.tiff"
writeRaster(mean_july_mean_max, filename=paste0(output_path, output_file), overwrite=TRUE, NAflag=-9999)






# plot
# 1) define breaks & colours
brks <- c(-Inf, 0.4, 0.5, 0.6, 1.8, 2.4, 3.1, Inf)
cols <- c("#fde725", "#dce319", "#b8de29", "#95d840", "#73d055", "#55c667", "#006400")

# 2) plot with classification
plot(mean_july_max,
     main = "Mean of Yearly Maximum July Interception Loss\n(Germany, 2005–2013)",
     xlab = "Longitude", ylab = "Latitude",
     asp = 1)


#######

# convert it to a SpatRaster
r_tm <- rast(mean_july_max)

# 2) Define breaks, labels & colours
breaks <- c(-Inf, 0.4, 0.5, 0.6, 1.8, 2.4, 3.1, Inf)
labs   <- c("≤0.4","0.4–0.5","0.5–0.6","0.6–1.8","1.8–2.4","2.4–3.1",">3.1")
cols   <- c("#fde725","#dce319","#b8de29","#95d840","#73d055","#55c667","#006400")

# 3) Build classification matrix (from, to, becomes)
m <- cbind(
  from    = breaks[-length(breaks)],
  to      = breaks[-1],
  becomes = seq_along(labs)
)

# 4) Classify into integer codes 1:7
r_cat_july <- classify(r_tm, rcl = m, include.lowest = TRUE, right = TRUE)

# 5) Turn into a factor
r_cat_july <- as.factor(r_cat_july)

# 6) Attach the full RAT so levels 1–7 exist even if some are empty
labs <- c("≤0.4","0.4–0.5","0.5–0.6","0.6–1.8","1.8–2.4","2.4–3.1",">3.1")
levels(r_cat_july) <- data.frame(ID = 1:7, category = labs)

# 7) Plot all 7 classes with your palette
cols <- c("#fde725","#dce319","#b8de29","#95d840","#73d055","#55c667","#006400")

plot(
  r_cat_july,
  col      = cols,
  legend   = TRUE,
  main     = "Mean July Interception Loss (Germany, 2005–2013)",
  xlab     = "Longitude",
  ylab     = "Latitude"
)
# manual legend with all 7 labels & colours
legend(
  "bottomright",
  legend = labs,
  fill   = cols,
  title  = "Loss (mm)",
  bty    = "n"
)
