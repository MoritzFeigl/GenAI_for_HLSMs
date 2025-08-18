
library(magrittr)
library(ggplot2)
library(wesanderson)
library(patchwork)
main_path <- "C:/Users/morit/Dropbox/Projekte/GenAI_para_for_HLSMs/analysis_results/validation_results"

# Create validation environment ----------------------------------------------------------

plot_data <- read.csv(paste0(main_path, "/validation_plot_data.csv"))
plot_data <- plot_data[, -1]
# select best runs, mHM run3, GenAI_para run 3, GenAI_para+SP run 1
plot_data <- plot_data[!(plot_data$run == "mHM" & plot_data$run == 3), ]
plot_data <- plot_data[!(plot_data$run == "GenAI_para" & plot_data$run == 3), ]
plot_data <- plot_data[!(plot_data$run == "GenAI_para + SP" & plot_data$run == 1), ]
plot_data <- plot_data[plot_data$experiment %in% c("mHM", "GenAI_para", "GenAI_para + SP"), ]
# get LSTM data
lstm_data <- read.csv("C:/Users/morit/Dropbox/Projekte/GenAI_para_for_HLSMs/results/06_regional_lstm/ensemble_run/test_ensemble_metrics.csv")
plot_data <- rbind(plot_data, 
                   data.frame("experiment"= "LSTM", "run"=1, "Basin"=lstm_data$basin,
                           "variable"="NSE", "value"=lstm_data$NSE, "tf_type"="LSTM"))
plot_data <- rbind(plot_data,
                   data.frame("experiment"= "LSTM", "run"=1, "Basin"=lstm_data$basin,
                           "variable"="lNSE", "value"=lstm_data$lNSE, "tf_type"="LSTM"))
plot_data <- rbind(plot_data, 
                   data.frame("experiment"= "LSTM", "run"=1, "Basin"=lstm_data$basin,
                           "variable"="KGE", "value"=lstm_data$KGE, "tf_type"="LSTM"))
# plot_data <- rbind(plot_data, 
#                    data.frame("experiment"= "LSTM", "run"=1, "Basin"=lstm_data$basin,
#                            "variable"="SPAEF", "value"=NA, "tf_type"="LSTM"))

# factors
plot_data$variable <- factor(plot_data$variable, 
                             levels=c("NSE", "lNSE", "KGE", "SPAEF"), 
                             labels=c("NSE", "log NSE", "KGE", "SPAEF"))
plot_data$experiment <- factor(plot_data$experiment, 
                               levels =c("mHM", "LSTM", "GenAI_para", "GenAI_para + SP"),
                               labels =c("Default-TFs\n(mHM)", 
                                         "LSTM\n(Kratzert et al. 2019)", 
                                         "GenAI-TFs v1\n(mHM)", 
                                         "GenAI-TFs\n(mHM)")) 
plot_data$tf_type <- factor(plot_data$tf_type, 
                            levels=c("mHM", "LSTM", "GenAI_para"), 
                            labels=c("mHM", "LSTM", "GenAI_para-mHM")) #


# p <- ggplot(plot_data[!(plot_data$variable %in% c("SPAEF", "log NSE")), ], 
#        aes(experiment, value)) + 
#   geom_boxplot() + 
#   # geom_jitter(width = 0.2, alpha=0.4) + 
#   facet_wrap(~variable, scales = "free", nrow=2)  + 
#   coord_cartesian(ylim = c(-1,1)) +
#   labs(x="", y="") +
#   theme_classic(base_family="sans", base_size=7) +
#   theme(
#     axis.title = element_text(size = 7),
#     axis.text = element_text(size = 7),
#     strip.text = element_text(size = 7),
#     strip.background = element_blank()
#     )

p1 <- ggplot(plot_data[plot_data$variable == "NSE" &
                         plot_data$experiment != "GenAI-TFs v1\n(mHM)", ], aes(experiment, value)) + 
  geom_boxplot() + 
  coord_cartesian(ylim = c(-2,1)) +
  labs(x="", y="Nash–Sutcliffe efficiency (NSE)") +
  theme_classic(base_family="sans", base_size=7) +
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    strip.text = element_text(size = 7),
    strip.background = element_blank()
  )
p1_dens <- ggplot(plot_data[plot_data$variable == "NSE" &
                         plot_data$experiment != "GenAI-TFs v1\n(mHM)", ], 
             aes(x = value, color = experiment)) + 
  stat_ecdf(geom = "step") + 
  coord_cartesian(xlim = c(-2, 1)) +
  labs(x = "Nash–Sutcliffe efficiency (NSE)", y = "Cumulative Probability") +
  theme_classic(base_family = "sans", base_size = 7) +
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    strip.text = element_text(size = 7),
    strip.background = element_blank(),
    legend.title = element_blank(),
    axis.title.x = element_text(vjust = 6)
  )

patchwork = p1 + p1_dens
patchwork = patchwork + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 12))


ggsave(
  filename = file.path(main_path, "validation_results_LSTM_NSE_cum.pdf"),
  plot     = patchwork,
  width    = 18,#8.8,
  height   = 7.5,#15,
  units    = "cm",
  dpi      = 300
)


p2 <- ggplot(plot_data[plot_data$variable == "KGE" &
                         plot_data$experiment != "GenAI-TFs v1\n(mHM)", ], aes(experiment, value)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(-2,1)) +
  labs(x="", y="Kling-Gupta efficiency (KGE)") +
  theme_classic(base_family="sans", base_size=7) +
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    strip.text = element_text(size = 7),
    strip.background = element_blank()
  )
p2_dens <- ggplot(plot_data[plot_data$variable == "KGE" &
                              plot_data$experiment != "GenAI-TFs v1\n(mHM)", ], 
                  aes(x = value, color = experiment)) + 
  stat_ecdf(geom = "step") + 
  coord_cartesian(xlim = c(-2, 1)) +
  labs(x = "Nash–Sutcliffe efficiency (NSE)", y = "Cumulative Probability") +
  theme_classic(base_family = "sans", base_size = 7) +
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    strip.text = element_text(size = 7),
    strip.background = element_blank(),
    legend.title = element_blank(),
    axis.title.x = element_text(vjust = 6)
  )
patchwork2 = p2 + p2_dens
patchwork2 = patchwork2 + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 12))


ggsave(
  filename = file.path(main_path, "validation_results_LSTM_KGE_cum.pdf"),
  plot     = patchwork2,
  width    = 18,
  height   = 7.5,
  units    = "cm",
  dpi      = 300
)


# 
# patchwork = p1 + p2
# patchwork = patchwork + plot_annotation(tag_levels = "a") & 
#   theme(plot.tag = element_text(face = "bold", size = 12))
# 
# 
# ggsave(
#   filename = file.path(main_path, "validation_results_LSTM.pdf"),
#   plot     = patchwork,
#   width    = 18,#8.8,    
#   height   = 7.5,#15,      
#   units    = "cm",
#   dpi      = 300      
# )
# ggsave(
#   filename = file.path(main_path, "validation_results_LSTM_NSE.pdf"),
#   plot     = p1,
#   width    = 8.8,
#   height   = 7.5,
#   units    = "cm",
#   dpi      = 300
# )
# ggsave(
#   filename = file.path(main_path, "validation_results_LSTM_KGE.pdf"),
#   plot     = p2,
#   width    = 8.8,
#   height   = 7.5,
#   units    = "cm",
#   dpi      = 300
# )



# get median values
aggregate(value ~ experiment + variable, plot_data, median)


