# This is for creating the table comparing the performances metrics for all models #
library(patchwork)
library(ggplot2)

# source all files in the R directory
files_to_source <- list.files("R/", recursive = TRUE, full.names = TRUE)
invisible(lapply(files_to_source, function(x) source(x, chdir = TRUE)))


# config setting
data_path_play_dyn_results <- "playground/2_Model_Structure_for_Baseline_and_Dynamic_Model/performances/DYN/"



init_DYN_results_plot <- function(){
  
  summary_plot <- tibble(metric = character(),
                         mean = numeric(),
                         lower.ci = numeric(),
                         upper.ci = numeric(),
                         median = numeric(),
                         min = numeric(),
                         max = numeric(), 
                         Q1 = numeric(),
                         Q3 = numeric(),
                         model = character(),
                         LM = numeric(),
                         imputation = character())
  
  return(summary_plot)
}
summary_plot <- init_DYN_results_plot()


datasets_files <- list.files(data_path_play_dyn_results, 
                             recursive = TRUE, full.names = TRUE)

datasets_files <- datasets_files[!grepl("trash", datasets_files)]


for (f in datasets_files){
  
  print(f)
  start_time <- Sys.time()
  
  # load data
  load(f) # loads results data
  
  results[is.infinite(results$value),]$value <- NA
  
  performance_summary <- results  %>% 
    group_by(model, metric, LM) %>% 
    mutate(mean = mean(value, na.rm = TRUE),
           sd = sd(value, na.rm = TRUE),
           n = n(),
           min = min(value, na.rm = TRUE),
           max = max(value, na.rm = TRUE),
           median = median(value, na.rm = TRUE),
           Q1 = quantile(value, probs = 0.25, na.rm = TRUE),
           Q3 = quantile(value, probs = 0.75, na.rm = TRUE),
           `median (IQR), %` = sprintf("%.3f (%.3f - %.3f)", median, Q1, Q3)) %>% 
    ungroup() %>%
    select(c("metric", "mean", "sd", "n", "min", "max", "median", "Q1", "Q3", "LM","model_type", "horizon", "imputation")) %>%
    unique() %>%
    mutate(se = sd / sqrt(n),
           lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se,
           upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se,
           `median (IQR), %` = sprintf("%.3f (%.3f , %.3f)", median, Q1, Q3),
           `mean (95% CI)` = sprintf("%.3f (%.3f , %.3f)", mean, lower.ci, upper.ci),
           model = paste(unique(model_type), unique(imputation), collapse = " "),
           LM = LM,
           imputation = imputation)
  
  
  summary_plot <- summary_plot %>%
    add_row(metric = performance_summary$metric,
            mean = performance_summary$mean,
            lower.ci = performance_summary$lower.ci,
            upper.ci = performance_summary$upper.ci,
            median = performance_summary$median,
            min = performance_summary$min,
            max = performance_summary$max,
            Q1 = performance_summary$Q1,
            Q3 = performance_summary$Q3,
            model = performance_summary$model,
            LM = performance_summary$LM,
            imputation = performance_summary$imputation) 
  
  summary_plot <- unique(summary_plot)
  
  
}

summary_plot <- summary_plot %>%
  mutate(model=if_else(model=="LM_COX MF", "LM Cox",
                       if_else(model=="LM_CS MF", "LM CS",
                               if_else(model=="LM_FG MF", "LM FG",
                                       if_else(model=="LM_FG_SEP MF", "FG-sep",
                                               if_else(model=="LM_LOG MF", "LM LR",
                                                       if_else(model=="LM_MULTI MF", "LM MLR",
                                                               if_else(model=="RMTL_TS MF", "RMTL-ts", model))))))))

plot1 <- summary_plot %>% 
  filter(metric == "AUROC") %>%
  mutate(model = as.factor(model)) %>%
  ggplot(aes(x = LM, y = mean, group = model, colour = model, linetype = model, shape = model)) +
  scale_x_continuous(breaks = c(0,7,14,21,28)) +
  geom_line(position = position_dodge(width = 0.1)) +
  scale_color_manual(values = c('#4A9E00', '#009E79','#F65246', '#AD8800', '#0091BD', '#8A66FF', '#FB41D0')) +
  scale_linetype_manual(values = c("solid", "dashed", "solid", "dotdash", "longdash", "twodash", "dashed")) +
  scale_shape_manual(values = c(1, 2, 3, 7, 5, 6, 4)) +
  geom_errorbar(aes(ymax = upper.ci, ymin = lower.ci), width = 0.1,
                position = position_dodge(width = 0.1)) +
  geom_point(position = position_dodge(width = 0.1)) +
  labs(y = "AUROC", title = "(a) AUC meanplot across different landmarks") + 
  theme_bw() +
  theme(
    # Hide panel borders and remove grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  coord_cartesian(ylim = c(0.593, 0.75))
  # + ylim(0.58, 0.75)


plot2 <- summary_plot %>% filter(metric == "slope") %>%
  mutate(model=as.factor(model)) %>%
  ggplot(aes(x = LM, y = mean, group = model, colour = model, linetype = model, shape = model)) +
  scale_x_continuous(breaks = c(0,7,14,21,28)) +
  geom_line(position = position_dodge(width = 0.1)) +
  scale_color_manual(values = c('#4A9E00', '#009E79','#F65246', '#AD8800', '#0091BD', '#8A66FF', '#FB41D0')) +
  scale_linetype_manual(values = c("solid", "dashed", "solid", "dotdash", "longdash", "twodash", "dashed")) +
  scale_shape_manual(values = c(1, 2, 3, 7, 5, 6, 4)) +
  geom_errorbar(aes(ymax = upper.ci, ymin = lower.ci), width = 0.1,
                position = position_dodge(width = 0.1)) +
  geom_point(position = position_dodge(width = 0.1)) +
  labs(y = "Calibration Slope", title = "(b) Calibration slope meanplot across different landmarks") + 
  theme_bw() +
  theme(
    # Hide panel borders and remove grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  coord_cartesian(ylim = c(0.348, 1.50)) +
  # ylim(0.348, 1.50) +
  geom_hline(yintercept = 1, linetype = "solid", color = "black") +
  annotate("text", x = Inf, y = 1.05, 
           label="ideal calibration slope = 1", hjust = 1.05, vjust = 0.5, size=3.5)


plot3 <- summary_plot %>% filter(metric == "OE_ratio") %>%
  mutate(model=as.factor(model)) %>%
  ggplot(aes(x = LM, y = mean, group = model, colour = model, linetype = model, shape = model)) +
  scale_x_continuous(breaks = c(0,7,14,21,28)) +
  geom_line(position = position_dodge(width = 0.1)) +
  scale_color_manual(values = c('#4A9E00', '#009E79','#F65246', '#AD8800', '#0091BD', '#8A66FF', '#FB41D0')) +
  scale_linetype_manual(values = c("solid", "dashed", "solid", "dotdash", "longdash", "twodash", "dashed")) +
  scale_shape_manual(values = c(1, 2, 3, 7, 5, 6, 4)) +
  geom_errorbar(aes(ymax = upper.ci, ymin = lower.ci), width = 0.1,
                position = position_dodge(width = 0.1)) +
  geom_point(position = position_dodge(width = 0.1)) +
  labs(y = "OE_ratio", title = "(c) O/E ratio meanplot across different landmarks") + 
  theme_bw() +
  theme(
    # Hide panel borders and remove grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  coord_cartesian(ylim = c(0.5, 1.4)) +
  # ylim(0.5, 1.4) +
  geom_hline(yintercept = 1, linetype = "solid", color = "black") +
  annotate("text", x = 12, y = 0.9, 
           label="ideal O/E ratio = 1", hjust = 1.5, vjust = 0.5, size=3.5)


plot4 <- summary_plot %>% filter(metric == "ECI") %>%
  mutate(model=as.factor(model)) %>%
  ggplot(aes(x = LM, y = mean, group = model, colour = model, linetype = model, shape = model)) +
  scale_x_continuous(breaks = c(0,7,14,21,28)) +
  geom_line(position = position_dodge(width = 0.1)) +
  scale_color_manual(values = c('#4A9E00', '#009E79','#F65246', '#AD8800', '#0091BD', '#8A66FF', '#FB41D0')) +
  scale_linetype_manual(values = c("solid", "dashed", "solid", "dotdash", "longdash", "twodash", "dashed")) +
  scale_shape_manual(values = c(1, 2, 3, 7, 5, 6, 4)) +
  geom_errorbar(aes(ymax = upper.ci, ymin = lower.ci), width = 0.1,
                position = position_dodge(width = 0.1)) +
  geom_point(position = position_dodge(width = 0.1)) +
  labs(y = "ECI", title = "(d) ECI meanplot across different landmarks") + 
  theme_bw() +
  theme(
    # Hide panel borders and remove grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  coord_cartesian(ylim = c(0.0, 0.19)) 
  # + ylim(0.0, 0.19)


plot5 <- summary_plot %>% filter(metric == "Scaled_BS") %>%
  mutate(model=as.factor(model)) %>%
  ggplot(aes(x = LM, y = mean, group = model, colour = model, linetype = model, shape = model)) +
  scale_x_continuous(breaks = c(0,7,14,21,28)) +
  geom_line(position = position_dodge(width = 0.1)) +
  scale_color_manual(values = c('#4A9E00', '#009E79','#F65246', '#AD8800', '#0091BD', '#8A66FF', '#FB41D0')) +
  scale_linetype_manual(values = c("solid", "dashed", "solid", "dotdash", "longdash", "twodash", "dashed")) +
  scale_shape_manual(values = c(1, 2, 3, 7, 5, 6, 4)) +
  geom_errorbar(aes(ymax = upper.ci, ymin = lower.ci), width = 0.1,
                position = position_dodge(width = 0.1)) +
  geom_point(position = position_dodge(width = 0.1)) +
  labs(y = "Scaled Brier Score", title = "(e) Scaled BS meanplot across different landmarks") + 
  theme_bw() +
  theme(
    # Hide panel borders and remove grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  coord_cartesian(ylim = c(-0.046, 0.03)) 
  # + ylim(-0.046, 0.03)

combined_plot <- plot1 + plot2 + plot3 + plot4 + plot5 + plot_layout(ncol = 2)

File <- paste0("playground/2_Model_Structure_for_Baseline_and_Dynamic_Model/Supplementary file 9/", "Fig3", ".tif")
tiff(File, width = 12, height = 14, units = "in", res = 300)

print(combined_plot)

dev.off()

