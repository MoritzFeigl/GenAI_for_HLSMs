setwd("C:/Users/morit/Dropbox/Projekte/GenAI_para_for_HLSMs/analysis_results/validation_results")
library(dplyr)

extract_run <- function(data, exp_name, run, var_prefix){
  sub <- data[data$experiment == exp_name & data$run == run, c("Basin", "variable", "value")]
  all_list <- list()
  for (var in c("NSE", "KGE", "lNSE")){
    tmp <- sub[sub$variable == var, c("Basin", "value")]
    names(tmp)[2] <- paste0(var_prefix, "-", var)
    all_list[[var]] <- tmp
  }
  extracted <- merge(merge(all_list[["NSE"]], all_list[["lNSE"]], by="Basin"), 
                     all_list[["KGE"]], by="Basin")
  return(extracted)
}

# Tables for GIS -------------------------------------------------------------------------
# Validation 
# mHM resultate
# Exp2-1 resultate
# 2, 1, 0.5 km resultate
km2 <- read.csv("validation_2km_plot_data.csv")
km1 <- read.csv("validation_1km_plot_data.csv")
km05 <- read.csv("validation_0.5km_plot_data.csv")


# GenAI_para Ergebnisse
GenAI_para <- read.csv("exp2-run1_val_results.csv", sep=";")[, 1:4]
names(GenAI_para) <- c("Basin", "GenAI_para-KGE", "GenAI_para-NSE", "GenAI_para-lNSE")

GenAI_para_2km <- extract_run(km2, "GenAI_para + SP", 1, "GenAI_para-2km")
GenAI_para_1km <- extract_run(km1, "GenAI_para + SP", 1, "GenAI_para-1km")
GenAI_para_05km <- extract_run(km05, "GenAI_para + SP", 1, "GenAI_para-0.5km")


# mHM Ergebnisse
mhm <- read.csv("exp0-run3_val_results.csv")[, 1:4]
names(mhm) <- c("Basin", "mHM-KGE", "mHM-NSE", "mHM-lNSE")
mhm_2km <- extract_run(km2, "mHM", 3, "mHM-2km")
mhm_1km <- extract_run(km1, "mHM", 3, "mHM-1km")
mhm_05km <- extract_run(km05, "mHM", 3, "mHM-0.5km")


merged_df <- GenAI_para %>%
  inner_join(GenAI_para_2km, by = "Basin") %>%
  inner_join(GenAI_para_1km, by = "Basin") %>%
  inner_join(GenAI_para_05km, by = "Basin") %>%
  inner_join(mhm, by = "Basin") %>%
  inner_join(mhm_2km, by = "Basin") %>%
  inner_join(mhm_1km, by = "Basin") %>%
  inner_join(mhm_05km, by = "Basin") 
merged_df$split <- "Validation"


# Training
training_results <- read.csv("../training_results_2024-01-16.csv")
GenAI_para_train <- training_results[training_results$experiment == 2 & training_results$run == 1,]
GenAI_para_train <- GenAI_para_train[, c("Basin", "NSE", "lNSE", "KGE")]
names(GenAI_para_train) <- c("Basin", "GenAI_para-KGE", "GenAI_para-NSE", "GenAI_para-lNSE")
mhm_train <- training_results[training_results$experiment == 0 & training_results$run == 3,]
mhm_train <- mhm_train[, c("Basin", "NSE", "lNSE", "KGE")]
names(mhm_train) <- c("Basin", "mHM-KGE", "mHM-NSE", "mHM-lNSE")

train_res <- merge(GenAI_para_train, mhm_train, by = "Basin")
for (var in names(merged_df)[c(5:13, 17:25)]){
  train_res <- cbind(train_res, new_var = NA)
  names(train_res)[length(names(train_res))] <- var
}
train_res$split <- "Training"
train_res <- train_res[, names(merged_df)]


all_results <- rbind(merged_df, train_res)


basin_properties <- read.table("../../data/01_study_basins_info/LUT_german_basins_with_properties_and_water_balance.txt", 
                               header=TRUE, sep=";")
col_sub <- c(1:7, 18:26, 29:34)
all <- merge(all_results, basin_properties[col_sub], by.x="Basin", by.y = "Stat_ID", 
             all.x = TRUE, all.y=FALSE)

all[which(is.na(all), arr.ind = TRUE)] <- -999

write.csv(all, "../compiled_gis_results.csv")



val <- all[all$split == "Validation", ]
summary(val$`GenAI_para-NSE` - val$`mHM-NSE`)

# Time series ----------------------------------------------------------------------------
setwd("C:/Users/morit/Dropbox/Projekte/GenAI_para_for_HLSMs/results/Validation/exp2-run1/")

basin <- "9316284"
ts <- read.table(paste0(basin, "_daily_discharge.out"), header = TRUE)
library(ggplot2)
ts$time = as.POSIXct(paste0(ts$Year, "-", ts$Mon, "-", ts$Day))
ts <- ts[, 5:7]
names(ts) <- c("Qobs", "Qsim", "time")
ts_melt <- reshape2::melt(ts, id.vars = "time")
ggplot(ts_melt[format(ts_melt$time, "%Y") %in% c("2011", "2012"), ], aes(x=time, y=value, col=variable)) + geom_line()






