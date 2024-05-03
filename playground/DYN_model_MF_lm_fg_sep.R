# this is for separate landmark fine-gray model #

# source all files in the R directory
files_to_source <- list.files("R/", recursive = TRUE, full.names = TRUE)
invisible(lapply(files_to_source, function(x) source(x, chdir = TRUE)))

# config for this model
model <- "DYN_MF_LM_FG_SEP"
model_type <- "LM_FG_SEP"
horizon <- "7 days"
imputation <- "MF"
path_data_complete <- data_path_play_dyn_complete_MF
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

# for time saving purpose, select only the first 10 list files
# train_files <- train_files[1:10]
# test_files <- test_files[1:10]

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
  
  # subset both train and test datasets with LM <= 30
  # around 5%-6% data were thrown away
  
  data_train_bin <- data_train_bin %>% filter(LM <= 30)
  data_test_bin <- data_test_bin %>% filter(LM <= 30)
  
  # apply administrative censoring to our dataset due to the nature of sliding prediction window
  # assume the prediction time horizon is set as 7 days
  
  # outcome: Surv(eventtime, type)
  
  data_train_bin <- data_train_bin %>%
    mutate(eventtime = ifelse(eventtime <= LM + 7, eventtime, LM + 7),
           type = ifelse(eventtime < LM + 7, type, "Censored"),
           Tstart = LM,
           id = paste(functioneelDossierNr, CAT_catheter_episode, sep = "_")) %>%
    filter(!eventtime <= Tstart) %>% 
    mutate(type = if_else(type == "Censored", 0, if_else(type == "CLABSI", 1, if_else(type == "Death", 2, 3))))
  
  data_test_bin <- data_test_bin %>%
    mutate(eventtime = ifelse(eventtime <= LM + 7, eventtime, LM + 7),
           type = ifelse(eventtime < LM + 7, type, "Censored"),
           Tstart = LM,
           id = paste(functioneelDossierNr, CAT_catheter_episode, sep = "_")) %>%
    filter(!eventtime <= Tstart) %>% 
    mutate(type = if_else(type == "Censored", 0, if_else(type == "CLABSI", 1, if_else(type == "Death", 2, 3))))
  
  
  # train data set
  
  Knots <- Hmisc::rcspline.eval(data_train_bin$CARE_VS_systolic_BP_last, nk=3, knots.only = TRUE) 
  CARE_VS_systolic_BP_last_nodes <- Hmisc::rcspline.eval(data_train_bin$CARE_VS_systolic_BP_last, nk=3) 
  attr(CARE_VS_systolic_BP_last_nodes, "dim") <- NULL
  attr(CARE_VS_systolic_BP_last_nodes, "knots") <- NULL
  
  data_train_bin$LAB_WBC_count_last_log <- log(data_train_bin$LAB_WBC_count_last)
  data_train_bin$LAB_CRP_last_log <- log(data_train_bin$LAB_CRP_last)
  
  data_train_bin$CARE_VS_systolic_BP_last_nodes <- CARE_VS_systolic_BP_last_nodes
  
  # test data set
  
  CARE_VS_systolic_BP_last_nodes <- Hmisc::rcspline.eval(data_test_bin$CARE_VS_systolic_BP_last, knots = c(Knots)) 
  attr(CARE_VS_systolic_BP_last_nodes, "dim") <- NULL
  attr(CARE_VS_systolic_BP_last_nodes, "knots") <- NULL
  data_test_bin$CARE_VS_systolic_BP_last_nodes <- CARE_VS_systolic_BP_last_nodes
  
  data_test_bin$LAB_WBC_count_last_log <- log(data_test_bin$LAB_WBC_count_last)
  data_test_bin$LAB_CRP_last_log <- log(data_test_bin$LAB_CRP_last)
  
  
  # specify the formula
  predictors_col <- vars_selected_2_BASE_DYN
  vars_not_in_model <- c("CARE_VS_systolic_BP_last", "LAB_WBC_count_last", "LAB_CRP_last")
  form <- as.formula(                      # Create formula
    paste(" ~ LAB_WBC_count_last_log + LAB_CRP_last_log + CARE_VS_systolic_BP_last_nodes +", paste0(predictors_col[!predictors_col %in% vars_not_in_model], collapse = " + ")))
  
  ############## not only coefficients, but also prediction metrics shall be saved separately #####
  library(prodlim)
  
  coefs_list <- list()
  results_list <- list()
  for(i in 0:30){
    data_train_i <- subset(data_train_bin, LM == i)
    fg_i <- riskRegression::FGR(update.formula(prodlim::Hist(eventtime,type)~., form), data_train_i, cause = 1)
    coefs_list[[paste("LM", i)]] <- fg_i$crrFit$coef
    data_test_i <- subset(data_test_bin, LM == i)
    data_test_i$pred_test <- predict(fg_i, newdata=data_test_i, times = i+7)
    results_list[[paste("LM", i)]] <- data_test_i[, c("pred_test", "eventtime", "type", "functioneelDossierNr", "CAT_catheter_episode")]
  }
  
  combined_results <-data.table::rbindlist(results_list, use.names = TRUE, fill = TRUE, idcol=TRUE)
  colnames(combined_results) <- c("LM", "pred_test", "eventtime", "type", "functioneelDossierNr", "CAT_catheter_episode")
  combined_results$y_test <- if_else(combined_results$type == 1, 1, 0)
  combined_results$LM <- substr(combined_results$LM, 3, 5)
  combined_results$LM <- as.numeric(combined_results$LM)
  combined_results$y_true_cat <- if_else(combined_results$type == 1, "CLABSI", if_else(combined_results$type == 2, "Death", if_else(combined_results$type == 3, "Discharge", "Censored")))
  
  
  coefs_results <- as.data.frame(do.call(rbind, coefs_list))
  coefs_results <- rownames_to_column(coefs_results, var = "LM")
  coefs_results$LM <- substr(coefs_results$LM, 3, 5)
  coefs_results$LM <- as.numeric(coefs_results$LM)
  coefs_results <- coefs_results %>% pivot_longer(cols= LAB_WBC_count_last_log:MS_is_ICU_unit,
                                                  names_to="variable",
                                                  values_to="estimates")
  
  
  predictions <- predictions %>% 
    add_row(preds = combined_results$pred_test,
            y_true_cat = combined_results$y_true_cat,
            y_true_time = combined_results$eventtime,
            train_set = str_replace(f, path_data_complete, ""),
            test_set = str_replace(test_file, path_data_complete, ""),
            model = model,
            model_type = model_type,
            horizon = horizon,
            imputation = imputation,
            LM = combined_results$LM,
            functioneelDossierNr = combined_results$functioneelDossierNr,
            CAT_catheter_episode = combined_results$CAT_catheter_episode)
  

 

  # performance metrics of the model
  
  metrics <- combined_results %>% group_by(LM) %>% summarise(
    AUROC = c(c_statistic(y_test, pred_test)),
    slope = calibration_slope(y_test, pred_test),
    OE_ratio = oe_ratio(y_test, pred_test),
    ECI = ECI(y_test, pred_test),
    Scaled_BS = scaled_brier_score(y_test, pred_test)) %>%
    pivot_longer(cols= AUROC:Scaled_BS,
                 names_to="metric",
                 values_to="value")
  
  results <- results %>%
    add_row(train_set = str_replace(f, path_data_complete, ""),
            test_set = str_replace(test_file, path_data_complete, ""),
            metric = metrics$metric,
            value = metrics$value,
            model = model,
            model_type = model_type,
            horizon = horizon,
            imputation = imputation,
            LM = metrics$LM)
  
  coefs <- coefs %>% 
    add_row(variable = coefs_results$variable,
            value = coefs_results$estimates,
            train_set = str_replace(f, path_data_complete, ""),
            model = model,
            model_type = model_type,
            horizon = horizon,
            imputation = imputation,
            LM = coefs_results$LM)
  
   message(sprintf("DONE in %s minutes.", 
                  difftime(Sys.time(), start_time, units = "mins") %>% as.numeric()))
}


# save predictions 
save(predictions, file = paste0("playground/2_Model_Structure_for_Baseline_and_Dynamic_Model/predictions/DYN/", preds_name))
save(results, file = paste0("playground/2_Model_Structure_for_Baseline_and_Dynamic_Model/performances/DYN/", results_name))
save(coefs, file = paste0("playground/2_Model_Structure_for_Baseline_and_Dynamic_Model/coefs/DYN/", coefs_name))
