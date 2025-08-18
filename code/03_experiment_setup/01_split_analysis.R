# GenAI_para-ET-mHM project
# Training/validation split for basins

if(Sys.info()[["nodename"]] == "cfgrammar"){
  setwd("/mnt/Data/Dropbox/Projekte/GenAI_para_for_HLSMs")
} else {
  if(Sys.info()["sysname"] == "Linux"){
    setwd("/home/lv71468/mfeigl/GenAI_para_for_HLSMs")
  } else {
    setwd("C:/Users/morit/Dropbox/Projekte/GenAI_para_for_HLSMs")
  }  
}
# libraries
library(ggplot2)


# 1. Get data ----------------------------------------------------------------------------
# Check available time series for all basins
study_basins <- read.csv("data/04_GIS/study_basins.csv")
study_basins <- sapply(study_basins$layer, function(x) unlist(strsplit(x, "_"))[3])
basins <- data.frame(Stat_ID = as.integer(study_basins),
                     Start_Date = NA, End_Date = NA, 
                     Easting_X=NA, Norting_Y=NA, River=NA, area = NA, dem_mean=NA, 
                     AET = NA, PET = NA, P = NA, Qobs = NA)
basin_properties <- read.table("data/01_study_basins_info/LUT_german_basins_with_properties_and_water_balance.txt", 
                               header=TRUE, sep=";")

for(id in seq_along(study_basins)){
  try({
    info <- readLines(paste0("data/03_discharge/", study_basins[id], ".txt"))
    rel_line <- grep("start", info)
    info <- info[rel_line:(rel_line+1)]
    basins[id, c("Start_Date", "End_Date")] <- c(substring(info[1], 8, 17),
                                                 substring(info[2], 8, 17))  
  })
  basins[id, c(4:12)] <- basin_properties[basin_properties$Stat_ID == basins$Stat_ID[id], 
                                          c("Easting_X", "Norting_Y", "River", 
                                            "catArea", "dem_mean", "AET.mm_a.1.", 
                                            "PET.mm_a.1.", "Precip.mm_a.1.", "Qobs.mm_a.1.")]
}


basins$Start_Date <- as.Date(basins$Start_Date, format="%Y %m %d")
basins$End_Date <- as.Date(basins$End_Date, format="%Y %m %d")
basins$start_year <- as.integer(format(basins$Start_Date, "%Y"))
basins$start_month <- as.integer(format(basins$Start_Date, "%m"))
basins$end_year <- as.integer(format(basins$End_Date, "%Y"))
basins$end_month <- as.integer(format(basins$End_Date, "%m"))

# add main river info
main_rivers <- read.csv("data/01_study_basins_info/main_river_info.csv")

basins <- merge(basins, main_rivers, all.x = TRUE)
# add names
basins$main_river_name <- ""
basins$main_river_name[basins$main_river == 6342800] <- "Donau"
basins$main_river_name[basins$main_river == 6340600] <- "Mulde"
basins$main_river_name[basins$main_river == 6335304] <- "Main"
basins$main_river_name[basins$main_river == 6335600] <- "Neckar"
basins$main_river_name[basins$main_river == 6337200] <- "Weser"
basins$main_river_name[basins$main_river == 6338100] <- "Ems"
basins$main_river_name[basins$main_river == 6340300] <- "Saale"
basins$main_river <- factor(basins$main_river)

# add main river names to the main river basins
main_river_id <- c(6342800, 6340600, 6335304, 6335600, 6337200, 6338100, 6340300)
main_river_name <- c("Donau", "Mulde", "Main", "Neckar", "Weser", "Ems", "Saale")
for(mriver in 1:7){
  basins$main_river[basins$Stat_ID == main_river_id[mriver]] <- main_river_id[mriver]
  basins$main_river_name[basins$Stat_ID == main_river_id[mriver]] <- main_river_name[mriver]
}

# 2. properties of training data ---------------------------------------------------------
# training 2014-2019

# training basins from random sampled basins with area < 1000 km²
sample <- read.csv("data/04_GIS/training_basins.csv")
sample <- sapply(sample$layer, function(x) as.integer(unlist(strsplit(x, "_"))[3]))

basins$split <- "Validation"
basins$split[basins$Stat_ID %in% sample] <- "Training"
basins$split[basins$end_year < 2008] <- "invalid"

# save split overview
dir.create("data/06_training_validation_data", showWarnings = FALSE)
write.csv(basins, "data/06_training_validation_data/basin_split.csv", row.names = FALSE)

# boxplots
basins_climate_data <- basins[,c("area", "dem_mean", "AET", "PET", "P", "Qobs", "split")]
basins_climate_data$log_area <- log(basins_climate_data$area)
basins_climate_data$area <- NULL
plot_data <- reshape2::melt(basins_climate_data, id.var = "split")
ggplot(plot_data,
       aes(x = split, y = value, fill = split)) + 
  geom_boxplot() + 
  geom_point(position = position_jitterdodge(), alpha = 0.4) +
  facet_wrap(. ~ factor(variable), scales="free_y") + 
  labs(x="", y="", fill="")
dir.create("results/01_split_analysis", showWarnings = FALSE)
ggsave("results/01_split_analysis/split_catchment_properties.png", width = 12, height = 8)

write.csv(basins, "results/01_split_analysis/basin_split.csv", row.names=FALSE)
