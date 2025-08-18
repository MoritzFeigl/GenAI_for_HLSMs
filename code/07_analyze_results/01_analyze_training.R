# GenAI_para-ET-mHM project
# analyze training

setwd("/gpfs/data/fs71468/GenAI_para_runs/results")
#setwd("C:/Users/morit/Dropbox/Projekte/GenAI_para_for_HLSMs/results")
dir.create("../analysis_results/training_results", recursive = TRUE)
library(magrittr)
library(ggplot2)
library(wesanderson)

# 1. all runs QCs ------------------------------------------------------------------------

training_results <- NULL
all_tfs <- NULL
for(experiment in c(0, 1, 2)){
  # for each experiment:
  #   1. for each run: plot time-series of best TF losses
  #   2. get overall best train results for all 50 training basins
  #   3. plot median (over 50 training basins) training results
  exp_result <- NULL
  for(run in 1:5){
    cat(experiment, "-", run, "\n")
    path <- paste0("exp", experiment, "-run", run, "/")
    # find current best result
    tf1 <- tf2 <- NULL
    try({tf1 <- read.csv(paste0(path, "best_TFs_tracker_1.csv"))})
    try({tf2 <- read.csv(paste0(path, "best_TFs_tracker_2.csv"))})
    last_tfs <- rbind(tf1[nrow(tf1), ], tf2[nrow(tf2), ])
    best_result <- last_tfs[which.max(c(tail(tf1$loss, 1), tail(tf2$loss, 1))), ]
    
    if(experiment != 0){
      best_tfs <- cbind(data.frame("experiment" = experiment, "run" = run),
            best_result[, 4:11])
      all_tfs <- rbind(all_tfs, best_tfs)
    }
    
    
    # 1. plot ts of best losses
    # plot some losses and best losses time series
    read_losses <- try({all_losses <- rbind(read.csv(paste0(path, "losses_tracker_1.csv")), 
                    read.csv(paste0(path, "losses_tracker_2.csv")))})
    if(class(read_losses) == "try-error"){
      loss_1 <- NULL
      loss_2 <- NULL
      try({loss_1 <- read.csv(paste0(path, "losses_tracker_1.csv"))})
      try({loss_2 <- read.csv(paste0(path, "losses_tracker_2.csv"))})
      all_losses <- rbind(loss_1, loss_2)
    }
    
    all_losses <- all_losses[!is.na(all_losses$loss), c("time", "loss")]
    all_losses$time <- as.POSIXct(all_losses$time)
    losses <- rbind(tf1, tf2)[, c("time", "loss")]
    losses$time <- as.POSIXct(losses$time)
    names(losses)[2] <- "best_loss"
    plot_losses <- merge(all_losses, losses, all=TRUE)
    current_best <- NA
    start_losses <- which(!is.na(plot_losses$best_loss))[1]
    current_best <- plot_losses[start_losses, "best_loss"]
    for(i in start_losses:nrow(plot_losses)){
      if(is.na(plot_losses[i, "loss"])){
        plot_losses[i, "best_loss"] <- current_best
        next
      }
      if(current_best > plot_losses[i, "loss"]){
        plot_losses[i, "best_loss"] <- current_best
      } else {
        current_best<- plot_losses[i, "loss"]
        plot_losses[i, "best_loss"] <- current_best
      }
    }
    
    plot_losses_melted <- reshape2::melt(plot_losses, id="time")
    plot_losses_melted$variable <- factor(plot_losses_melted$variable, levels=c("best_loss", "loss"))
    
    experiment_labels <- c("0" = "mHM", "1" = "GenAI_para", "2" = "GenAI_para + SP")
    ggplot(plot_losses_melted, aes(time, value, col=variable)) +
      geom_line()+
      scale_color_manual(name="", labels=c('Best model loss', 'Iteration Loss'), 
          values=c("black", alpha("grey", 0.7))) +
      coord_cartesian(ylim=c(-1, 1)) +
      labs(x="", y="Loss", title=paste0(experiment_labels[as.character(experiment)], " - run ", run))
    ggsave(paste0("../analysis_results/training_results/", 
    "exp", experiment, "-run", run, "_optimization_progression.png"),
           width = 14, height = 8, units = "in")
    
    # get best run results
    qc1 <- qc2 <- NULL
    qc1 <- try({feather::read_feather(paste0(path, "basin_qc_tracker_1.feather"))})
    qc2 <- try({feather::read_feather(paste0(path, "basin_qc_tracker_2.feather"))})
    qc <- rbind(qc1, qc2)
    
    best_qc <- qc[qc$time == best_result$time, ]
    
    exp_result <- cbind(data.frame("experiment" = experiment, "run" = run,
                                   "mean_loss" = best_result$loss),
                        best_qc)
    training_results <- rbind(training_results, exp_result)
  }
}
table(training_results$run)
date <- Sys.Date()
write.csv(all_tfs, paste0("best_tfs_", date, ".csv"), row.names = FALSE)
write.csv(training_results, paste0("training_results_", date, ".csv"), row.names = FALSE)

# Boxplots -------------------------------------------------------------------------------

training_results <- read.csv(paste0("training_results_", date, ".csv"))

# 2 versions: all basins, without bad basins
bad_data <- training_results[training_results$NSE < 0, ]
bad_basins <- names(table(bad_data$Basin))[table(bad_data$Basin) >25]

# 1. Plot each quality criteria for each experiment and run as boxplots
plot_data <- reshape2::melt(training_results[, c("experiment", "Basin", "run", "KGE", "NSE", "lNSE", "SPAEF")],
                            id.vars = c("experiment", "Basin", "run"))
plot_data$experiment <- factor(plot_data$experiment, levels = c(0, 1, 2),
                               labels = c("mHM", "GenAI_para", "GenAI_para + SP"))
for(variable in c("NSE", "lNSE", "KGE", "SPAEF")){
  ggplot(plot_data[plot_data$variable == variable, ], 
         aes(factor(experiment), value, fill=factor(run))) + 
    geom_boxplot() + 
    scale_fill_manual(values=wes_palette(name="Darjeeling1")) +
    labs(x="Experiment", fill="Run", y = variable)
  ggsave(paste0("../analysis_results/training_results/all_training_runs_", variable, ".png"),
         width = 10, height = 8, units = "in")
}


plot_data <- plot_data[!(plot_data$Basin %in% bad_basins), ]
for(variable in c("NSE", "lNSE", "KGE", "SPAEF")){
  ggplot(plot_data[plot_data$variable == variable, ], 
         aes(factor(experiment), value, fill=factor(run))) + 
    geom_boxplot() + 
    scale_fill_manual(values=wes_palette(name="Darjeeling1")) +
    labs(x="Experiment", fill="Run", y = variable)
  ggsave(paste0("../analysis_results/training_results/all_training_runs_", variable, "_NoBadBasins.png"),
         width = 10, height = 8, units = "in")
}

# 2. Plot only best runs
best_losses <- aggregate(mean_loss ~ experiment, training_results, max)
best_runs <- unique(merge(best_losses, training_results[, c("experiment", "run", "mean_loss")], 
                          all.x = TRUE, all.y = FALSE))
best_training <- merge(best_runs, training_results, all.x=TRUE, all.y=FALSE)

# function to get median values in boxplots
get_box_stats <- function(y, upper_limit=1) {
  return(data.frame(
    y = upper_limit,
    label = paste(
      "Median =", round(median(y), 2)
    )
  ))
}


plot_data <- reshape2::melt(best_training[, c("experiment", "Basin", "run", "KGE", "NSE", "lNSE", "SPAEF")],
                            id.vars = c("experiment", "Basin", "run"))
plot_data$tf_type <- "GenAI_para"
plot_data$tf_type[plot_data$experiment == 0] <- "mHM"
plot_data$tf_type <- factor(plot_data$tf_type, levels = c("GenAI_para", "mHM"))

plot_data$experiment <- factor(plot_data$experiment, levels = c(0, 1, 2),
                               labels = c("mHM", "GenAI_para", "GenAI_para + SP"))
plot_data$variable <- factor(plot_data$variable, levels=c("NSE", "lNSE", "KGE", "SPAEF"))

# with SPAEF
ggplot(plot_data, 
       aes(factor(experiment), value, fill=tf_type)) + 
  geom_boxplot() + 
  geom_jitter(width = 0.2, alpha=0.4) + 
  facet_wrap(~variable, scales = "free") + 
  scale_fill_manual("TF type", values=wes_palette(name="Darjeeling1", n = 4)[c(2,4)]) +
  labs(x="Experiment", fill="TF type", y = "") +
  stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0, size=3)
ggsave(paste0(main_path, "analysis_results/training_results/best_training_runs.png"),
       width = 18, height = 8, units = "in")


# without spaef
ggplot(plot_data[plot_data$variable != "SPAEF", ], 
       aes(factor(experiment), value, fill=tf_type)) + 
  geom_boxplot() + 
  geom_jitter(width = 0.2, alpha=0.4) + 
  facet_wrap(~variable, scales = "free") + 
  scale_fill_manual("TF type", values=wes_palette(name="Darjeeling1", n = 4)[c(2,4)]) +
  labs(x="Experiment", fill="TF type", y = "") +
  stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0, size=3)
ggsave("../analysis_results/training_results/best_training_runs_NoSpaef.png",
       width = 18, height = 8, units = "in")



















