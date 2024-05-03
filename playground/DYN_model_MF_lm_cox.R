# source all files in the R directory
files_to_source <- list.files("R/", recursive = TRUE, full.names = TRUE)
invisible(lapply(files_to_source, function(x) source(x, chdir = TRUE)))

# config for this model
model <- "DYN_MF_LM_COX"
model_type <- "LM_COX"
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
           type = ifelse(eventtime < LM + 7 & type == "CLABSI", 1, 0),
           Tstart = LM,
           id = paste(functioneelDossierNr, CAT_catheter_episode, sep = "_")) %>%
    filter(!eventtime <= Tstart) 
  
  data_test_bin <- data_test_bin %>%
    mutate(eventtime = ifelse(eventtime <= LM + 7, eventtime, LM + 7),
           type = ifelse(eventtime < LM + 7 & type == "CLABSI", 1, 0),
           Tstart = LM,
           id = paste(functioneelDossierNr, CAT_catheter_episode, sep = "_")) %>%
    filter(!eventtime <= Tstart) 
  

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
  
  
  # also create s, s^2, ICU*s, ICU*s^2
  data_train_bin$LM_1 <- data_train_bin$LM/30
  data_train_bin$LM_2 <- (data_train_bin$LM/30)^2
  data_train_bin$MS_is_ICU_unit_1 <-  data_train_bin$MS_is_ICU_unit * data_train_bin$LM_1
  data_train_bin$MS_is_ICU_unit_2 <-  data_train_bin$MS_is_ICU_unit * data_train_bin$LM_2
  

  # also create s, s^2, ICU*s, ICU*s^2
  data_test_bin$LM_1 <- data_test_bin$LM/30
  data_test_bin$LM_2 <- (data_test_bin$LM/30)^2
  data_test_bin$MS_is_ICU_unit_1 <-  data_test_bin$MS_is_ICU_unit * data_test_bin$LM_1
  data_test_bin$MS_is_ICU_unit_2 <-  data_test_bin$MS_is_ICU_unit * data_test_bin$LM_2
  
  predictors_col <- c("CAT_catheter_type_binary_all_CVC", 
                      "CAT_catheter_type_binary_all_Port_a_cath",
                      "CAT_catheter_type_binary_all_Tunneled_CVC", 
                      "CAT_catheter_type_binary_all_PICC", 
                      "CAT_catheter_location_binary_all_Collarbone",
                      "CAT_catheter_location_binary_all_Neck",
                      "MED_7d_TPN", 
                      "MED_L2_7d_J01_ANTIBACTERIALS_FOR_SYSTEMIC_USE", 
                      "MED_L2_7d_L01_ANTINEOPLASTIC_AGENTS",
                      "CLABSI_history",
                      "COM_PATH_tumor_before_LM",
                      "COM_lymphoma_before_LM",
                      "COM_PATH_transplant_before_LM",
                      "MB_other_infection_than_BSI_during_window",
                      "ADM_admission_source_binary_all_Home",
                      "CARE_VS_MV",
                      "MS_is_ICU_unit",
                      "CARE_VS_temperature_max",
                      "CARE_VS_systolic_BP_last_nodes",
                      "LAB_WBC_count_last_log",
                      "LAB_CRP_last_log",
                      "LM_1",
                      "LM_2",
                      "MS_is_ICU_unit_1",
                      "MS_is_ICU_unit_2")
  
  # specify the formula
  form <- paste0(" ~ cluster(id) + ",paste0(predictors_col, collapse = " + "))
  

  # fit landmark supermodel
  # simple
  # use LMsupercox0
  LMsupercox0 <- coxph(update.formula(Surv(Tstart,eventtime,type==1)~., form), data_train_bin, method="breslow", x=TRUE)
  
  coef_cox_model <- coef(LMsupercox0)
  
  # predict on test set
 
  bh <- basehaz(LMsupercox0,centered=FALSE)
  lp <- predict(LMsupercox0, newdata = data_test_bin, reference = "zero")
  
  data_test_bin$pred_test <- NA
  
  tt <- data_test_bin$LM
  nt <- length(data_test_bin$LM)

  bh$haz0 <- c(0, diff(bh$hazard))
  bh$Haz0 <- cumsum(bh$haz0)
  
  landmark=7
  Lambda0<-bh$Haz0
  etimes<-bh$time
  
  # in case vector is too big
  p <- do.call("cbind",lapply(tt, function(tti){
    index <- prodlim::sindex(eval.times=c(tti,landmark+tti),jump.times=etimes)
    Lambda0.diff <- c(0,Lambda0)[1+index[2]] - c(0,Lambda0)[1+index[1]]
  }))
  
  data_test_bin$pred_test=1-exp(-as.vector(p) * exp(lp))

  
  # observed risk within 7 days y_true_cat (categorical: "CLABSI", "no_CLABSI", "Discharge", "Death", "Censored")
  # observed risk within 7 days y_test (binary: 0/1)
  data_test_bin$y_true_cat <- if_else(data_test_bin$type == 1, "CLABSI", "Censored")
  
  data_test_bin$y_test <- data_test_bin$type
  
  
  predictions <- predictions %>% 
    add_row(preds = data_test_bin$pred_test,
            y_true_cat = data_test_bin$y_true_cat,
            y_true_time = data_test_bin$eventtime,
            train_set = str_replace(f, path_data_complete, ""),
            test_set = str_replace(test_file, path_data_complete, ""),
            model = model,
            model_type = model_type,
            horizon = horizon,
            imputation = imputation,
            LM = data_test_bin$LM,
            functioneelDossierNr = data_test_bin$functioneelDossierNr,
            CAT_catheter_episode = data_test_bin$CAT_catheter_episode)
  
  # performance metrics of the model
  
  metrics <- data_test_bin %>% group_by(LM) %>% summarise(
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
    add_row(variable = names(coef_cox_model),
            value = coef_cox_model,
            train_set = str_replace(f, path_data_complete, ""),
            model = model,
            model_type = model_type,
            horizon = horizon,
            imputation = imputation,
            LM = NA)
  
  message(sprintf("DONE in %s minutes.", 
                  difftime(Sys.time(), start_time, units = "mins") %>% as.numeric()))
}

# save predictions 
save(predictions, file = paste0("playground/2_Model_Structure_for_Baseline_and_Dynamic_Model/predictions/DYN/", preds_name))
save(results, file = paste0("playground/2_Model_Structure_for_Baseline_and_Dynamic_Model/performances/DYN/", results_name))
save(coefs, file = paste0("playground/2_Model_Structure_for_Baseline_and_Dynamic_Model/coefs/DYN/", coefs_name))
