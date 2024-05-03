# source all files in the R directory
files_to_source <- list.files("R/", recursive = TRUE, full.names = TRUE)
invisible(lapply(files_to_source, function(x) source(x, chdir = TRUE)))

# config for this model
model <- "BASE_MF_FG_ALL_DAYS"
model_type <- "fg"
horizon <- "all days"
imputation <- "MF"
path_data_complete <- data_path_play_base_complete_MF
# selected variables
vars_selected_2_BASE_DYN <- vars_selected_DOMK_LIT


preds_name <- paste("preds", model, model_type, horizon, imputation, sep = "_")
results_name <- paste("results", model, model_type, horizon, imputation, sep = "_")
coefs_name <- paste("coefs", model, model_type, horizon, imputation, sep = "_")

# get filenames for imputed datasets 
datasets_files <- list.files(path_data_complete, 
                             recursive = TRUE, full.names = TRUE)
train_files <- datasets_files[str_detect(datasets_files, "train")]
test_files <- datasets_files[str_detect(datasets_files, "test")]

# keep results
predictions <- init_preds_BASE_DYN()
results <- init_results_BASE_DYN()
coefs <- init_coefs()

# build model for each df
for (f in train_files){
  
  print(f)
  start_time <- Sys.time()
  
  # load data
  load(f) # loads train data named data_train
  test_file <- str_replace(f, "train", "test") # corresponding test set file
  load(test_file) 
  
  # binarize
  cat_cols <- cat_cols_base
  
  bin_model <- data_train %>% 
    make_col_binary_drop(cat_cols)
  data_train_bin <- bin_model$data
  data_test_bin <- data_test %>% 
    make_col_binary_drop(cat_cols, dropped_levels = bin_model$dropped_levels)
  data_test_bin <- data_test_bin$data
  
  # outcome: Surv(eventtime, type)
  # Surv(time, event)
  # time: follow-up time for right-censored data
  # event: status indicator 0=right censored, 1=event at time
  data_train_bin <-data_train_bin 
  
  data_test_bin <- data_test_bin 
  
  predictors_col <- vars_selected_2_BASE_DYN
  
  # specify the formula
  vars_not_in_model <- c("CARE_VS_systolic_BP_last", "LAB_WBC_count_last", "LAB_CRP_last")
  # fix the position of knots for CARE_VS_systolic_BP_last
  Knots <- Hmisc::rcspline.eval(data_train_bin$CARE_VS_systolic_BP_last, nk=3, knots.only = TRUE) 
  CARE_VS_systolic_BP_last_nodes <- Hmisc::rcspline.eval(data_train_bin$CARE_VS_systolic_BP_last, nk=3) 
  attr(CARE_VS_systolic_BP_last_nodes, "dim") <- NULL
  attr(CARE_VS_systolic_BP_last_nodes, "knots") <- NULL
  data_train_bin$CARE_VS_systolic_BP_last_nodes <- CARE_VS_systolic_BP_last_nodes
  
  data_train_bin$LAB_WBC_count_last_log <- log(data_train_bin$LAB_WBC_count_last)
  data_train_bin$LAB_CRP_last_log <- log(data_train_bin$LAB_CRP_last)
  
  
  # specify the formula
  form <- as.formula(                      # Create formula
    paste(" ~ LAB_WBC_count_last_log + LAB_CRP_last_log + CARE_VS_systolic_BP_last_nodes +", paste0(predictors_col[!predictors_col %in% vars_not_in_model], collapse = " + ")))
  # build model
  fg_model <- riskRegression::FGR(update.formula(prodlim::Hist(eventtime,type)~., form), data_train_bin, cause = "CLABSI")
  coef_fg_model <- fg_model$crrFit$coef
  
  # predict on test set
  library(prodlim)
  
  CARE_VS_systolic_BP_last_nodes <- Hmisc::rcspline.eval(data_test_bin$CARE_VS_systolic_BP_last, knots = c(Knots)) 
  attr(CARE_VS_systolic_BP_last_nodes, "dim") <- NULL
  attr(CARE_VS_systolic_BP_last_nodes, "knots") <- NULL
  data_test_bin$CARE_VS_systolic_BP_last_nodes <- CARE_VS_systolic_BP_last_nodes
  
  data_test_bin$LAB_WBC_count_last_log <- log(data_test_bin$LAB_WBC_count_last)
  data_test_bin$LAB_CRP_last_log <- log(data_test_bin$LAB_CRP_last)
  # predicted risk: pred_test
  
  data_test_bin$pred_test <- NA
  data_test_bin$pred_test <- predict(fg_model, newdata=data_test_bin, times = 7)


  # observed risk within 7 days y_true_cat (categorical: "CLABSI", "no_CLABSI", "Discharge", "Death", "Censored")
  # observed risk within 7 days y_test (binary: 0/1)
  data_test_bin$y_true_cat <-  if_else(data_test_bin$type == "CLABSI" & data_test_bin$eventtime <= 7, "CLABSI", 
                                       if_else(data_test_bin$type == "Death" & data_test_bin$eventtime <= 7, "Death",
                                               if_else(data_test_bin$type == "Discharge" & data_test_bin$eventtime <= 7, "Discharge", "Censored")))       
  data_test_bin$y_true_time <-  if_else(data_test_bin$type == "CLABSI" & data_test_bin$eventtime <= 7, data_test_bin$eventtime, 
                                       if_else(data_test_bin$type == "Death" & data_test_bin$eventtime <= 7, data_test_bin$eventtime,
                                               if_else(data_test_bin$type == "Discharge" & data_test_bin$eventtime <= 7, data_test_bin$eventtime, 7)))       
  
  data_test_bin$y_test <- if_else(data_test_bin$type == "CLABSI" & data_test_bin$eventtime <= 7, 1, 0)
  
  
  predictions <- predictions %>% 
    add_row(preds = data_test_bin$pred_test,
            y_true_cat = data_test_bin$y_true_cat,
            y_true_time = data_test_bin$y_true_time,
            train_set = str_replace(f, path_data_complete, ""),
            test_set = str_replace(test_file, path_data_complete, ""),
            model = model,
            model_type = model_type,
            horizon = horizon,
            imputation = imputation,
            LM = 0,
            functioneelDossierNr = data_test_bin$functioneelDossierNr,
            CAT_catheter_episode = data_test_bin$CAT_catheter_episode)
  
  
  # performance metrics of the model
  
  
  metrics <- c("AUROC" = c_statistic(data_test_bin$y_test, data_test_bin$pred_test))
  metrics["slope"] <- calibration_slope(data_test_bin$y_test, data_test_bin$pred_test)
  metrics["O/E"] <- oe_ratio(data_test_bin$y_test, data_test_bin$pred_test)
  metrics["ECI"] <- ECI(data_test_bin$y_test, data_test_bin$pred_test)
  metrics["Scaled_BS"] <- scaled_brier_score(data_test_bin$y_test, data_test_bin$pred_test)
  
  
  
  results <- results %>%
    add_row(train_set = str_replace(f, path_data_complete, ""),
            test_set = str_replace(test_file, path_data_complete, ""),
            metric = names(metrics),
            value = metrics,
            model = model,
            model_type = model_type,
            horizon = horizon,
            imputation = imputation,
            LM = 0)
  
  coefs <- coefs %>% 
    add_row(variable = names(coef_fg_model),
            value = coef_fg_model,
            train_set = str_replace(f, path_data_complete, ""),
            model = model,
            model_type = model_type,
            horizon = horizon,
            imputation = imputation,
            LM = 0)
  
  message(sprintf("DONE in %s minutes.", 
                  difftime(Sys.time(), start_time, units = "mins") %>% as.numeric()))
}

# save predictions 
save(predictions, file = paste0("playground/2_Model_Structure_for_Baseline_and_Dynamic_Model/predictions/BASE/", preds_name))
save(results, file = paste0("playground/2_Model_Structure_for_Baseline_and_Dynamic_Model/performances/BASE/", results_name))
save(coefs, file = paste0("playground/2_Model_Structure_for_Baseline_and_Dynamic_Model/coefs/BASE/", coefs_name))
##############################################################################################################################
load("C:/CLABSI/Garbage/clabsi_r/playground/2_Model_Structure_for_Baseline_and_Dynamic_Model/predictions/BASE/preds_BASE_MF_FG_ALL_DAYS_fg_all days_MF")
load("C:/CLABSI/Garbage/clabsi_r/playground/2_Model_Structure_for_Baseline_and_Dynamic_Model/performances/BASE/results_BASE_MF_FG_ALL_DAYS_fg_all days_MF")

library(pROC)
library(colorspace)
library(cutpointr)

# make AUROC curves
predictions_sub <- predictions %>% 
  mutate(y_test = ifelse(y_true_cat == "CLABSI", 1, 0)) %>%
  filter(!(is.na(preds)|preds>=1))

# summarize the AUROCs across 100 datasets
AUROC_summary <- results  %>% 
  filter(metric == "AUROC") %>% 
  summarise(mean_AUROC = mean(value, na.rm = TRUE),
            sd_AUROC = sd(value, na.rm = TRUE),
            n_AUROC = n(),
            median_AUROC = median(value),
            Q1_AUROC = quantile(value, probs = 0.25),
            Q3_AUROC = quantile(value, probs = 0.75)) %>%
  mutate(se_AUROC = sd_AUROC / sqrt(n_AUROC),
         lower.ci.AUROC = mean_AUROC - qt(1 - (0.05 / 2), n_AUROC - 1) * se_AUROC,
         upper.ci.AUROC = mean_AUROC + qt(1 - (0.05 / 2), n_AUROC - 1) * se_AUROC,
         `median (IQR), %` = sprintf("%.3f (%.3f - %.3f)", median_AUROC, Q1_AUROC, Q3_AUROC),
         `mean (95% CI)` = sprintf("%.3f (%.3f - %.3f)", mean_AUROC, lower.ci.AUROC, upper.ci.AUROC))


# add mean ROC using cutpointr for specifying the thresholds manually via the oc_manual function
# apply the same sequence of thresholds to all samples and take the mean of the sensitivity and specificity per threshold to get the "mean ROC curve"
# other option: summary ROC curves (SROC) for fitting a parametric model which combines multiple ROC curves

mr <- mean_roc(predictions_sub)

# ggplot(mr, aes(x = 1 - specificity, y = sensitivity)) + 
#   geom_step() + geom_point() +
#   theme(aspect.ratio = 1)
File <- paste0("playground/2_Model_Structure_for_Baseline_and_Dynamic_Model/plots/", preds_name, ".jpeg")
jpeg(File, quality = 100, width = 1400, height = 1400, units = "px")

cutpointr(data = predictions_sub, 
          x = preds, class = y_test, subgroup = test_set,
          pos_class = 1, direction = ">=") %>% 
  plot_roc(display_cutpoint = F) + theme(legend.position="none") +
  geom_line(data = mr, mapping = aes(x = 1 - specificity, y = sensitivity), 
            color = "black", linewidth = 1) +
  ggtitle(paste0(preds_name)) +
  theme(text = element_text(size = 24)) +
  theme(axis.text = element_text(size = 24)) +
  annotate(geom = "label", x = 0.8, y = 0.2, size = 10, label = paste(paste0("median (IQR), %: ", AUROC_summary$`median (IQR), %`), paste0("mean (95% CI): ", AUROC_summary$`mean (95% CI)`), sep="\n")) 

dev.off()
# make calibration plots

# summarize the calibration slopes and intercepts across 100 datasets

Calibrate_summary <- results  %>% 
  group_by(model, metric) %>% 
  mutate(mean = mean(value, na.rm = TRUE),
         sd = sd(value, na.rm = TRUE),
         n = n(),
         median = median(value),
         Q1 = quantile(value, probs = 0.25),
         Q3 = quantile(value, probs = 0.75),
         `median (IQR), %` = sprintf("%.3f (%.3f - %.3f)", median, Q1, Q3)) %>% 
  ungroup() %>%
  select(c("metric", "mean", "sd", "n", "median", "Q1", "Q3")) %>%
  unique() %>%
  filter(metric == "intercept"|metric == "slope") %>%
  mutate(se = sd / sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se,
         upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se,
         `median (IQR), %` = sprintf("%.3f (%.3f - %.3f)", median, Q1, Q3),
         `mean (95% CI)` = sprintf("%.3f (%.3f - %.3f)", mean, lower.ci, upper.ci))


Pal = qualitative_hcl(100)

data_list <- split(predictions_sub, f = predictions_sub$train_set) 

File <- paste0("playground/2_Model_Structure_for_Baseline_and_Dynamic_Model/plots/", results_name, ".jpeg")
jpeg(File, quality = 100, width = 1400, height = 1400, units = "px")

par(mar= c(5, 5, 4, 2) + 0.1)
plot(NULL,xlim=c(0,1),ylim=c(0,1),
     xlab="Predicted probability",ylab="Actual probability", cex.lab=2, cex.axis=2)
title(results_name, adj = 0, cex.main = 2.5, font.main = 1)

for(i in 1:100){
  
  predictions = data_list[[i]]$preds
  labels = data_list[[i]]$y_test
  logit_p = qlogis(predictions)
  logit_p[is.infinite(logit_p)] <- NA
  mod = rms:::lrm(labels~logit_p)
  pp = seq(0.01, 0.99, 0.01)
  pred = plogis(predict(mod, newdata = list(logit_p = qlogis(pp))))
  lines(pp, pred, col=Pal[i])
}

pred_overall = predictions_sub$preds
label_overall = predictions_sub$y_test
logit_p = qlogis(pred_overall)
logit_p[is.infinite(logit_p)] <- NA
mod = rms:::lrm(label_overall~logit_p)
pp = seq(0.01, 0.99, 0.01)
pred = plogis(predict(mod, newdata = list(logit_p = qlogis(pp))))
lines(pp, pred, col="black", lwd = 3)
legend(0.05, 0.8, col = 1:4, legend = c(paste0("mean slope (95% CI): ", Calibrate_summary$`mean (95% CI)`[1]), paste0("mean intercept (95% CI): ", Calibrate_summary$`mean (95% CI)`[2]), paste0("median slope (IQR), %: ", Calibrate_summary$`median (IQR), %`[1]), paste0("median intercept (IQR), %: ", Calibrate_summary$`median (IQR), %`[2])),
       xjust = 0.1,      # 0.5 means center adjusted
       yjust = 0.1,      # 0.5 means center adjusted
       x.intersp = -0.1, # adjust character interspacing as you like to effect box width
       y.intersp = 0.9,  # adjust character interspacing to effect box height
       adj = c(0.0, 0.0), cex = 2.5)

dev.off()

# this is the 10 group decile calibration plot
data_list <- split(predictions_sub, f = predictions_sub$train_set) 

File <- paste0("playground/2_Model_Structure_for_Baseline_and_Dynamic_Model/plots/", paste0(results_name, "02"), ".jpeg")
jpeg(File, quality = 100, width = 1400, height = 1400, units = "px")

par(mar= c(5, 5, 4, 2) + 0.1)
lim <- c(0, 0.1)
plot(NULL, xlim=lim, ylim=lim,
     xlab="Predicted probability",ylab="Actual probability", cex.lab=2, cex.axis=2)
title(results_name, adj = 0, cex.main = 2.5, font.main = 1)

# title(results_name, adj = 0, cex.main = 2.5, font.main = 1)

for(i in 1:100){
  
  predictions = data_list[[i]]$preds
  labels = data_list[[i]]$y_test
  table_s = cbind(as.data.frame(predictions), as.data.frame(labels))
  colnames(table_s) <- c("yhat", "y")
  table_s <- lift_table(table_s, bin_number = 10)
  lines(table_s$mean_yhat, table_s$mean_y, col="grey", pch = 17, cex=1.5, type = "b", lty = 2)
}

pred_overall = predictions_sub$preds
label_overall = predictions_sub$y_test
table_s = cbind(as.data.frame(pred_overall), as.data.frame(label_overall))
colnames(table_s) <- c("yhat", "y")
table_s <- lift_table(table_s, bin_number = 10)
lines(table_s$mean_yhat, table_s$mean_y, col="black", pch = 17, cex=3, type = "b", lty = 2, lwd = 3)
abline(coef = c(0,1), lwd = 2)

rcs_obj <- rcspline_plot(pred_overall, label_overall, model="logistic", nk=3, show="prob", plot = FALSE)
lines(rcs_obj$x, rcs_obj$xbeta, col="blue", lwd = 3)
# intervals
lines(rcs_obj$x, rcs_obj$upper, lty = 'dashed', col = 'blue', lwd = 3)
lines(rcs_obj$x, rcs_obj$lower, lty = 'dashed', col = 'blue', lwd = 3)


sm = lowess(pred_overall, label_overall, iter=0)
lines(sm, col='red', lwd = 3)

x <- pred_overall
bins <- seq(lim[1], lim[2], length=101)
x <- x[x >= lim[1] & x <= lim[2]]
f <- table(cut(x, bins))
j <- f > 0
bins <- (bins[-101])[j]
f <- f[j]
f <- lim[1] + .15 * diff(lim) * f / max(f)
segments(bins, 0, bins, f)

legend("topleft",  legend=c("Ideal", "Flexible calibration (RCS)", "Flexible calibration (Loess)", "Grouped observations (g = 10)"),
       col=c("black", "blue", "red", "black"), lty = c(1,1,1,2), lwd = c(2,3,3,3), xjust = 0.1, yjust = 0.1, cex = 2.5)

dev.off()

