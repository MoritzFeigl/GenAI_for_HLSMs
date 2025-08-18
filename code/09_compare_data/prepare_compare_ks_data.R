setwd("C:/Users/morit/Dropbox/Projekte/GenAI_para_for_HLSMs/data/10_compare_data/Hydraul_Param_SoilGrids_Schaap_0")

library(ncdf4)
library(raster)
library(sf)
library(rnaturalearth)

# 1) Read the fill value from the NetCDF
ncfile <- "Hydraul_Param_SoilGrids_Schaap_sl4.nc"
nc     <- nc_open(ncfile)
fillv  <- ncatt_get(nc, "mean_Ks_30cm", "_FillValue")$value
nc_close(nc)

# 2) Load the 15cm mean Ks raster and mask out no‐data
r <- raster(ncfile, varname = "mean_Ks_30cm")
r[r == fillv] <- NA

# 3) Download Germany outline via rnaturalearth
germany_sf <- ne_countries(country = "Germany", returnclass = "sf")

# 4) Crop & mask the raster to the Germany polygon
#    Convert sf → Spatial for use with raster::mask
germany_sp <- as(germany_sf, "Spatial")
r_de      <- crop(r, germany_sp)
r_de      <- mask(r_de, germany_sp)


# 5) Define breaks, labels and colors
breaks <- c(-Inf, 14, 25, 80, 115, 190, 220, Inf)
labs   <- c("≤14","14–25","25–80","80–115","115–190","190–220",">220")
cols   <- colorRampPalette(rev(c("brown","yellow","green","darkgreen")))(length(breaks)-1)
vals    <- getValues(r_de)  # raw numeric vector
fct_vals<- cut(
  vals,
  breaks        = breaks,
  include.lowest= TRUE,
  right         = TRUE,
  labels        = labs
)
# create an empty raster and assign integer codes
r_cat <- raster(r_de)
r_cat[] <- as.integer(fct_vals)

# attach a RAT (Raster Attribute Table) so plot.factor() knows the labels
rat <- data.frame(ID = 1:length(labs), category = labs)
levels(r_cat) <- rat

# 7) Plot
plot(
  r_cat,
  col    = cols,
  legend = TRUE,
  main   = "Mean Saturated Hydraulic Conductivity at 30 cm Depth (Germany)",
  xlab   = "Longitude",
  ylab   = "Latitude"
)



########### Dai et al. 2019
setwd("C:/Users/morit/Dropbox/Projekte/GenAI_para_for_HLSMs/data/10_compare_data/Dai_et_al_2019_Ks")

library(ncdf4)
library(raster)
library(sf)
library(rnaturalearth)

library(terra)

# 1) Open the NetCDF as a SpatRaster (no data are read yet)
r_all <- rast("k_s_l1.nc")

# 2) Define a simple lon/lat bbox for Germany
#    (xmin, xmax, ymin, ymax)
germ_ext <- ext(5, 15, 47, 55)

# 3) Crop—terra will compute which file-blocks intersect and only read those
r_de <- crop(r_all, germ_ext)

# 4) Mask to the precise border if you need
germany_sf <- ne_countries(country="Germany", returnclass="sf")
r_de <- mask(r_de, vect(germany_sf))

# 5) Replace fill and plot
r_de[r_de == -9999] <- NA


# 6) Define your breaks, labels & colors
breaks <- c(-Inf, 14, 25, 80, 115, 190, 220, Inf)
labs    <- c("≤14","14–25","25–80","80–115","115–190","190–220",">220")
cols    <- c("#c7552a","#f0daa5","#fbf2c4","#b8cdab", "#73a690", "#008383", "#004242")

# 7) Build a classification matrix: from, to, becomes
m <- matrix(c(
  -Inf,  14, 1,
  14,  25, 2,
  25,  80, 3,
  80, 115, 4,
  115, 190, 5,
  190, 220, 6,
  220, Inf, 7
), ncol=3, byrow=TRUE)

# 8) Classify into integer codes 1:7
r_cat <- classify(r_de, rcl = m, include.lowest = TRUE, right = TRUE)

# 9) Turn it into a factor and attach your labels
r_cat <- as.factor(r_cat)
levels(r_cat) <- data.frame(ID = 1:7, category = labs)

# 10) Plot with your custom palette
plot(
  r_cat,
  col      = cols,
  legend   = TRUE,
  main     = "Mean Saturated Hydraulic Conductivity at 30 cm Depth (Germany)",
  xlab     = "Longitude",
  ylab     = "Latitude"
)


########### Gupta et al. 2021
setwd("C:/Users/morit/Dropbox/Projekte/GenAI_para_for_HLSMs/data/10_compare_data/Gupta_et_al_2021")

library(ncdf4)
library(raster)
library(sf)
library(rnaturalearth)

library(terra)

# 1) Open the NetCDF as a SpatRaster (no data are read yet)
r_all <- rast("Global_Ksat_1Km_s0....0cm_v1.0.tif")

# 2) Define a simple lon/lat bbox for Germany
#    (xmin, xmax, ymin, ymax)
germ_ext <- ext(5, 15, 47, 55)

# 3) Crop—terra will compute which file-blocks intersect and only read those
r_de <- crop(r_all, germ_ext)

# 4) Mask to the precise border if you need
germany_sf <- ne_countries(country="Germany", returnclass="sf")
r_de <- mask(r_de, vect(germany_sf))

# 5) Replace fill and plot
r_de[r_de == -9999] <- NA

# 3) Invert the log10 transform
r_de <- 10 ^ r_de

# 8) Export raster to GeoTIFF
crs(r_de) <- "EPSG:4326"
output_path <- "C:/Users/morit/Dropbox/Projekte/GenAI_para_for_HLSMs/analysis_results/compare_tiffs/"
dir.create(output_path, showWarnings = FALSE)
output_file <- "Ksat_Germany_Gupta_et_al_20221.tiff"
writeRaster(r_de, filename=paste0(output_path, output_file), overwrite=TRUE, NAflag=-9999)


# 6) Define your breaks, labels & colors
breaks <- c(-Inf, 14, 25, 80, 115, 190, 220, Inf)
labs    <- c("≤14","14–25","25–80","80–115","115–190","190–220",">220")
cols    <- c("#c7552a","#f0daa5","#fbf2c4","#b8cdab", "#73a690", "#008383", "#004242")

# 7) Build a classification matrix: from, to, becomes
m <- matrix(c(
  -Inf,  14, 1,
  14,  25, 2,
  25,  80, 3,
  80, 115, 4,
  115, 190, 5,
  190, 220, 6,
  220, Inf, 7
), ncol=3, byrow=TRUE)

# 8) Classify into integer codes 1:7
r_cat <- classify(r_de, rcl = m, include.lowest = TRUE, right = TRUE)

# 9) Turn it into a factor and attach your labels
r_cat <- as.factor(r_cat)
levels(r_cat) <- data.frame(ID = 1:7, category = labs)

# 10) Plot with your custom palette
plot(
  r_cat,
  col      = cols,
  legend   = TRUE,
  main     = "Mean Saturated Hydraulic Conductivity at 30 cm Depth (Germany)",
  xlab     = "Longitude",
  ylab     = "Latitude"
)

