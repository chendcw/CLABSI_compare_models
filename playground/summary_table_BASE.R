# This is for creating the table comparing the performances metrics for all models #

# source all files in the R directory
files_to_source <- list.files("R/", recursive = TRUE, full.names = TRUE)
invisible(lapply(files_to_source, function(x) source(x, chdir = TRUE)))


# config setting
data_path_play_base_results <- "playground/2_Model_Structure_for_Baseline_and_Dynamic_Model/performances/BASE/"


summary_results <- init_results_sumstats()
summary_plot <- init_results_plot()


# get filenames for calculated predictions and results
datasets_files <- list.files(data_path_play_base_results, 
                             recursive = TRUE, full.names = TRUE)

datasets_files <- datasets_files[!grepl("trash", datasets_files)]


for (f in datasets_files){
  
  print(f)
  start_time <- Sys.time()
  
  # load data
  load(f) # loads results data
  
  results[is.infinite(results$value),]$value <- NA
  
  performance_summary <- results  %>% 
    group_by(model, metric) %>% 
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
    select(c("metric", "mean", "sd", "n", "min", "max", "median", "Q1", "Q3", "model_type", "horizon", "imputation")) %>%
    unique() %>%
    mutate(se = sd / sqrt(n),
           lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se,
           upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se,
           `median (IQR), %` = sprintf("%.3f (%.3f , %.3f)", median, Q1, Q3),
           `mean (95% CI)` = sprintf("%.3f (%.3f , %.3f)", mean, lower.ci, upper.ci),
           model = paste(unique(model_type), unique(horizon), collapse = " "),
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
            imputation = performance_summary$imputation)
  
  summary_plot <- unique(summary_plot)
  
  summary_results <- summary_results %>% 
    add_row(mean_AUC = performance_summary[performance_summary$metric=="AUROC",]$`mean (95% CI)`,
            mean_Calibration_Slope = performance_summary[performance_summary$metric=="slope",]$`mean (95% CI)`,
            mean_OE_ratio = performance_summary[performance_summary$metric=="O/E",]$`mean (95% CI)`,
            mean_ECI = performance_summary[performance_summary$metric=="ECI",]$`mean (95% CI)`,
            mean_Scaled_BS = performance_summary[performance_summary$metric=="Scaled_BS",]$`mean (95% CI)`,
            model = performance_summary$model,
            imputation = performance_summary$imputation)
  
  summary_results <- unique(summary_results)
  
}



# Write to table

library(xlsx)
write.xlsx(summary_results, file = "playground/2_Model_Structure_for_Baseline_and_Dynamic_Model/sumstatsBASE.xlsx")

