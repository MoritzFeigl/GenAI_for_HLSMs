# Plot paper results
library(ggplot2)
library(patchwork)
setwd("C:/Users/morit/Dropbox/Projekte/GenAI_para_for_HLSMs/analysis_results/validation_results/")
main_path <- "C:/Users/morit/Dropbox/Projekte/GenAI_para_for_HLSMs/analysis_results/validation_results"


# Best validation results
plot_data <- read.csv("validation_plot_data.csv")
experiments <- c("mHM", "GenAI_para", "GenAI_para + SP")
plot_data$experiment <- factor(plot_data$experiment, levels = experiments,
                               labels = experiments)
plot_data$variable <- factor(plot_data$variable, levels=c("NSE", "lNSE", "KGE", "SPAEF"))
write.csv(plot_data, paste0(main_path, "/analysis_results/validation_results/validation_0.5km_plot_data.csv"))





# Results on different sales

plot_data <- cbind(read.csv("validation_plot_data.csv"), "Resolution"="4x4 km")
plot_data2km <- cbind(read.csv("validation_2km_plot_data.csv"), "Resolution"="2x2 km")
plot_data1km <- cbind(read.csv("validation_1km_plot_data.csv"), "Resolution"="1x1 km")
plot_data05km <- cbind(read.csv("validation_0.5km_plot_data.csv"), "Resolution"="0.5x0.5 km")

spatial_plot <- rbind(cbind(read.csv("validation_plot_data.csv"), "Resolution"="4x4 km"),
                      cbind(read.csv("validation_2km_plot_data.csv"), "Resolution"="2x2 km"),
                      cbind(read.csv("validation_1km_plot_data.csv"), "Resolution"="1x1 km"),
                      cbind(read.csv("validation_0.5km_plot_data.csv"), "Resolution"="0.5x0.5 km"))



experiment_subset <- c("GenAI_para + SP")
spatial_plot$experiment <- factor(spatial_plot$experiment, levels = experiments,
                               labels = experiments)
spatial_plot$variable <- factor(spatial_plot$variable, 
                             levels=c("NSE", "lNSE", "KGE", "SPAEF"), 
                             labels=c("NSE", "log NSE", "KGE", "SPAEF"))


p1 <- ggplot(spatial_plot[spatial_plot$variable == "NSE" &
                           plot_data$experiment %in% experiment_subset, ], 
            aes(factor(Resolution), value, )) + 
  geom_boxplot() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x="Model Resolution", y="NSE") +
  theme_classic(base_family="sans", base_size=7) +
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    strip.text = element_text(size = 7),
    strip.background = element_blank()
  )
p2 <- ggplot(spatial_plot[spatial_plot$variable == "KGE" &
                            plot_data$experiment %in% experiment_subset, ], 
             aes(factor(Resolution), value, )) + 
  geom_boxplot() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x="Model Resolution", y="KGE") +
  theme_classic(base_family="sans", base_size=7) +
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    strip.text = element_text(size = 7),
    strip.background = element_blank()
  )
patchwork = p1 + p2
patchwork = patchwork + plot_annotation(tag_levels = "a") & 
                                        theme(plot.tag = element_text(face = "bold", size = 12))

ggsave(
  filename = file.path(main_path, "resolution_paper_results.pdf"),
  plot     = patchwork,
  width    = 18,    
  height   = 7.5,      
  units    = "cm",
  dpi      = 300      
)













