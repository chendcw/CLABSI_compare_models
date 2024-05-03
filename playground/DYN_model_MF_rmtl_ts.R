# source all files in the R directory
files_to_source <- list.files("R/", recursive = TRUE, full.names = TRUE)
invisible(lapply(files_to_source, function(x) source(x, chdir = TRUE)))

library(RMTL)

# config for this model
# config for this model; follows the naming convention: 
# [BASE/DYN]_[model type]_[var_sel]_[imputation]_[outcome type]_[horizon]_[other...]
# model name should be unique across files!!!
model <- "DYN_RMTL_LITDOMK_TS_MF_bin_7d"
print(model)
path_data_complete <- data_path_play_dyn_complete_MF

# model specific config
max_LM <- 30 
positive_class <- "CLABSI"
n_folds <- 5

# model specific functions
data_X_t <- function(data, t, center = TRUE, scale = TRUE){
  data %>% 
    filter(LM == t) %>% 
    select(-LM) %>% 
    scale(center = center, scale = scale) 
} 

data_Y_t <- function(data, t){
  data %>% filter(LM == t) %>% pull(CLABSI)
}

# cross-validation function modified
cvMTL_2 <- function(X, Y, type="Classification", Regularization="L21",
                    Lam1_seq=10^seq(1,-4, -1), Lam2=0, G=NULL, k=2,
                    opts=list(init=0, tol=10^-3, maxIter=1000),
                    stratify=FALSE, nfolds=5, ncores=2, parallel=FALSE){
  
  pretty <- function(y) pmax(pmin(y, 1 - 10^-15), 10^-15)
  
  calcError_2 <- function(m, newX=NULL, newY=NULL){
    if(class(m)!="MTL"){
      stop("The first arguement is not a MTL model")}
    if(!is.null(newX) & !is.null(newY)){
      task_num <- length(newY)
      yhat <- RMTL:::predict.MTL(m,newX)
      if(m$type=="Classification"){
        #residue <- lapply(1:task_num, function(x)
        #  newY[[x]]-(round(yhat[[x]])-0.5)*2)
        newY <- sapply(newY, function(x) ifelse(x == -1, 0, 1))
        # MSE
        # error <- sapply(1:task_num, function(x)
        #   mean((newY[[x]]-yhat[[x]])^2))
        # logloss
        error <- sapply(1:task_num, function(x) -mean(newY[[x]] * log(pretty(yhat[[x]])) + 
                                                        (1 - newY[[x]]) * log(1 - pretty(yhat[[x]]))))
        
      }else if(m$type=="Regression"){
        error <- sapply(1:task_num, function(x)
          mean((newY[[x]]-yhat[[x]])^2))
      }
      return(mean(error))
    }else{stop(" no new data (X or Y) are provided ")}
  }
  
  #test vilidity of input data
  if (!missing(X) & !missing(Y)){
    if (all(sapply(X, class)!="matrix")){
      X <- lapply(X, function(x){as.matrix(x)})
    }
    if (all(sapply(Y, class)!="matrix")){
      Y <- lapply(Y, function(x){as.matrix(x)})
    }
  }else{
    stop("data X or Y doesnot exists")
  }
  task_num <- length(X)
  if(stratify & type=="Regression"){
    stop("stratified CV is not applicable to regression")}
  cvPar <- RMTL:::getCVPartition(Y, nfolds, stratify)
  
  #cv
  if (!parallel){
    cvm <- rep(0, length(Lam1_seq))
    for (i in 1:nfolds){
      cv_Xtr <- lapply(c(1:task_num),
                       function(x) X[[x]][cvPar[[i]][[1]][[x]], ])
      cv_Ytr <- lapply(c(1:task_num),
                       function(x) Y[[x]][cvPar[[i]][[1]][[x]]])
      cv_Xte <- lapply(c(1:task_num),
                       function(x) X[[x]][cvPar[[i]][[2]][[x]], ])
      cv_Yte <- lapply(c(1:task_num),
                       function(x) Y[[x]][cvPar[[i]][[2]][[x]]])
      opt <- opts
      for (p_idx in 1: length(Lam1_seq)){
        m <- RMTL::MTL(X=cv_Xtr, Y=cv_Ytr, type=type,
                       Regularization=Regularization, Lam1=Lam1_seq[p_idx],
                       Lam2=Lam2, opts=opt, k=k, G=G)
        #non sparse model training
        if (!is.element(Regularization, c("Graph", "CMTL"))){
          opt$init <- 1
          opt$W0 <- m$W
          opt$C0 <- m$C
        }
        cv_err <- calcError_2(m, newX=cv_Xte, newY=cv_Yte)
        cvm[p_idx] = cvm[p_idx]+cv_err
      }
    }
    cvm = cvm/nfolds
  } else {
    requireNamespace('doParallel')
    requireNamespace('foreach')
    # modify to close cluster
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    #doParallel::registerDoParallel(ncores)
    cvm <- foreach::foreach(i = 1:nfolds, .combine="cbind") %dopar%{
      cv_Xtr <- lapply(c(1:task_num),
                       function(x) X[[x]][cvPar[[i]][[1]][[x]], ])
      cv_Ytr <- lapply(c(1:task_num),
                       function(x) Y[[x]][cvPar[[i]][[1]][[x]]])
      cv_Xte <- lapply(c(1:task_num),
                       function(x) X[[x]][cvPar[[i]][[2]][[x]], ])
      cv_Yte <- lapply(c(1:task_num),
                       function(x) Y[[x]][cvPar[[i]][[2]][[x]]])
      opt <- opts
      cvVec=rep(0, length(Lam1_seq))
      for (p_idx in 1: length(Lam1_seq)){
        m <- RMTL::MTL(X=cv_Xtr, Y=cv_Ytr, type=type,
                       Regularization=Regularization, Lam1=Lam1_seq[p_idx],
                       Lam2=Lam2, opts=opt, k=k, G=G)
        #non sparse model training
        if (!is.element(Regularization, c("Graph", "CMTL")) ){
          opt$init <- 1
          opt$W0 <- m$W
          opt$C0 <- m$C
        }
        cv_err <- calcError_2(m, newX=cv_Xte, newY=cv_Yte)
        cvVec[p_idx] <- cv_err
      }
      return(cvVec)
    }
    # modify to close cluster
    stopCluster(cl)
    cvm <- rowMeans(cvm)
  }
  
  best_idx <- which(cvm==min(cvm))[1]
  cv <- list(Lam1_seq=Lam1_seq, Lam1.min=Lam1_seq[best_idx],
             Lam2=Lam2, cvm=cvm)
  class(cv) <- "cvMTL"
  return(cv)
}


# build file names
preds_name <- paste("preds", model, sep = "_")
var_imp_name <- paste("var_imp", model, sep = "_")
hyperparams_name <- paste("hyperparams", model, sep = "_")
results_name <- paste("results", model, sep = "_")
timings_name <- paste("timings", model, sep = "_")

preds_file <- paste0(preds_dir_dyn, preds_name)
var_imp_file <- paste0(var_imp_dir_dyn, var_imp_name)
hyperparams_file <- paste0(models_path, hyperparams_name)
results_file <- paste0(preds_dir_dyn, results_name)
timings_file <- paste0(models_path, timings_name)

# get filenames for imputed datasets 
datasets_files <- list.files(path_data_complete, 
                             recursive = TRUE, full.names = TRUE)
train_files <- datasets_files[str_detect(datasets_files, "train")]
test_files <- datasets_files[str_detect(datasets_files, "test")]

# keep results
predictions <- init_preds()
var_imp <- init_var_imp()
hyperparams <- init_hyperparams()
timings <- init_timings()

OOB_predictions <- init_preds()

positive_class <- "CLABSI"
negative_class <- "no_CLABSI"

set.seed(2024)

# build model for each df
for (f in train_files){
  
  print(f)
  start_time <- Sys.time()
  
  predictions_this_model <- init_preds()
  
  # load data
  load(f) # loads train data named data_train
  test_file <- str_replace(f, "train", "test") # corresponding test set file
  load(test_file) 
  
  # filter LMs for which eventtime is at exact LM time
  # this is due to bug in the data prep - will be fixed for training 
  data_train <- data_train %>% 
    filter(eventtime > LM) 
  
  data_test <- data_test %>% 
    filter(eventtime > LM)
  
  # keep till LM 30
  data_train <- data_train %>% 
    filter(LM <= max_LM) 
  data_test <- data_test %>% 
    filter(LM <= max_LM) 
  
  # outcome 7 days, coded as -1, 1
  data_train <- data_train %>% 
    mutate(time_to_LM = eventtime - LM,
           CLABSI = if_else(type == "CLABSI" & time_to_LM <= 7
                            , 1, -1)) %>% 
    select(-c(eventtime, type, time_to_LM))
  
  data_test <- data_test %>% 
    mutate(time_to_LM = eventtime - LM,
           CLABSI = if_else(type == "CLABSI" & time_to_LM <= 7
                            , 1, -1)) %>% 
    select(-c(eventtime, type, time_to_LM))
  
  cols_keep <- vars_selected_DOMK_LIT[vars_selected_DOMK_LIT %in% colnames(data_train)]
  cols_keep <- cols_keep[!cols_keep %in% c("functioneelDossierNr", "type", 
                                           "eventtime")]
  cols_keep <- c("LM", cols_keep)
  
  # keep X and Y 
  train_X <- data_train %>% 
    select(all_of(cols_keep)) 
  
  test_X <- data_test %>% 
    select(all_of(cols_keep)) 
  
  # non linearities train data set
  Knots <- Hmisc::rcspline.eval(train_X$CARE_VS_systolic_BP_last, nk=3, knots.only = TRUE) 
  CARE_VS_systolic_BP_last_nodes <- Hmisc::rcspline.eval(train_X$CARE_VS_systolic_BP_last, nk=3) 
  attr(CARE_VS_systolic_BP_last_nodes, "dim") <- NULL
  attr(CARE_VS_systolic_BP_last_nodes, "knots") <- NULL
  
  train_X$LAB_WBC_count_last_log <- log(train_X$LAB_WBC_count_last)
  train_X$LAB_CRP_last_log <- log(train_X$LAB_CRP_last)
  
  train_X$CARE_VS_systolic_BP_last_nodes <- CARE_VS_systolic_BP_last_nodes
  
  train_X <- train_X %>% 
    select(-c(CARE_VS_systolic_BP_last, LAB_WBC_count_last, LAB_CRP_last))
  
  # non linearities test data set
  CARE_VS_systolic_BP_last_nodes <- Hmisc::rcspline.eval(test_X$CARE_VS_systolic_BP_last, knots = c(Knots)) 
  attr(CARE_VS_systolic_BP_last_nodes, "dim") <- NULL
  attr(CARE_VS_systolic_BP_last_nodes, "knots") <- NULL
  test_X$CARE_VS_systolic_BP_last_nodes <- CARE_VS_systolic_BP_last_nodes
  
  test_X$LAB_WBC_count_last_log <- log(test_X$LAB_WBC_count_last)
  test_X$LAB_CRP_last_log <- log(test_X$LAB_CRP_last)
  
  test_X <- test_X %>% 
    select(-c(CARE_VS_systolic_BP_last, LAB_WBC_count_last, LAB_CRP_last))
  
  # make sure columns are in the same order (for standardizing correctly)
  test_X <- test_X %>% 
    select(all_of(colnames(train_X)))
  
  Y_train <- data_train %>% 
    select(LM, CLABSI)
  Y_test <- data_test %>% 
    select(LM, CLABSI)
  
  # split per LM
  # train sets X
  data_X_train <- lapply(0:max_LM, function(t) data_X_t(train_X,t))
  # keep attribs
  scale_attribs <- lapply(data_X_train, attributes)
  scale_attribs_center <- lapply(scale_attribs, function(x) x[["scaled:center"]])
  scale_attribs_scale <- lapply(scale_attribs, function(x) x[["scaled:scale"]])
  # remove attribs
  data_X_train <- lapply(data_X_train, function(df) as.matrix(as.data.frame(df)))
  
  # train sets Y
  data_Y_train <- lapply(0:max_LM, function(t) data_Y_t(Y_train,t))
  
  # test sets X
  data_X_test <- lapply(0:max_LM, function(t) data_X_t(test_X,t, 
                                                       center = scale_attribs_center[[t+1]],
                                                       scale = scale_attribs_scale[[t+1]]))
  # test sets Y
  data_Y_test <- lapply(0:max_LM, function(t) data_Y_t(Y_test,t))
  
  # no. of events in train sets
  # lapply(data_Y_train, table)
  
  # check that standardizing is applied correctly
  # names(scale_attribs_center[[1]]) == colnames(data_X_test[[1]])
  
  # create graph
  G <- matrix(0, nrow = max_LM + 1, ncol = max_LM + 1)
  diag(G) <- 1
  G[row(G) - col(G) == 1] <- -1
  
  # build model
  start_time_tuning <- Sys.time()
  cv_fit <- cvMTL_2(data_X_train, data_Y_train, type="Classification",
                    Regularization="Graph",
                    Lam1_seq = 10^seq(1,-4, -0.1),
                    Lam2 = 0,
                    nfolds = n_folds,
                    G = G,
                    ncores = 12,
                    parallel = TRUE)
  time_tuning <- as.numeric(difftime(Sys.time(), start_time_tuning, units = "secs"))
  
  # train
  start_time_final_model <- Sys.time()
  model_MTL <- MTL(data_X_train, data_Y_train, type="Classification", 
                   Regularization="Graph",
                   G = G,
                   Lam1=cv_fit$Lam1.min,
                   Lam2 = 0)
  time_final_model <- as.numeric(difftime(Sys.time(), start_time_final_model, units = "secs"))
  
  # predict on test set
  start_time_predict <- Sys.time()
  test_preds <- predict(model_MTL, data_X_test)
  time_predict <- as.numeric(difftime(Sys.time(), start_time_predict, units = "secs"))
  
  test_preds_LM <- lapply(1:length(test_preds), function(i)
    tibble(value = test_preds[[i]][,1],
           LM = i)) %>% 
    do.call("rbind", .) %>% 
    mutate(LM = LM - 1)
  y_true_LM <- lapply(1:length(data_Y_test), function(i)
    tibble(value = data_Y_test[[i]],
           LM = i)) %>% 
    do.call("rbind", .) %>% 
    mutate(LM = LM - 1)
  
  predictions_this_model <- predictions_this_model %>% 
    add_row(preds = test_preds_LM$value,
            y_true_cat = ifelse(y_true_LM$value == 1, "CLABSI", "no_CLABSI"),
            y_true_time = NA_real_,
            train_set = str_replace(f, path_data_complete, ""),
            test_set = str_replace(test_file, path_data_complete, ""),
            model = model,
            LM = test_preds_LM$LM)
  
  # put back the functioneel dossier nr
  predictions_this_model <- predictions_this_model %>% 
    select(-c(functioneelDossierNr, CAT_catheter_episode))
  
  predictions_this_model <- data_test %>% 
    filter(LM <= max_LM) %>% arrange(LM) %>% 
    select(functioneelDossierNr, CAT_catheter_episode, CLABSI) %>% 
    cbind(predictions_this_model) %>% 
    as_tibble()
  
  # check that the stiching was correct
  n_mismatches <- predictions_this_model %>% 
    mutate(CLABSI = if_else(CLABSI == 1, "CLABSI", "no_CLABSI")) %>% 
    filter(CLABSI != y_true_cat) %>% 
    nrow()
  
  if(n_mismatches > 0){
    stop("Somthing went wrong when adding back functioneelDossierNr.")
  }
  
  predictions_this_model <- predictions_this_model %>% 
    select(-CLABSI)
  
  predictions <- predictions %>% 
    add_row(predictions_this_model)
  
  # model_weights
  model_weights <- model_MTL$W %>% 
    as.data.frame() %>% as_tibble() 
  
  cols_LMs <- colnames(model_weights)
  
  model_weights <- model_weights %>% 
    add_column(feat = colnames(data_X_train[[1]])) %>% 
    pivot_longer(cols = all_of(cols_LMs), names_to = "LM", values_to = "coef") %>% 
    mutate(LM = as.integer(str_replace(LM, "V", "")),
           LM = LM - 1) %>% # 1 is LM 0  
    group_by(feat) %>% 
    mutate(always_0 = all(coef == 0)) %>% 
    ungroup() %>% 
    filter(!always_0) %>% 
    select(-always_0)
  
  # these are standardized coeffs
  # to compare across LMs, we'd better unstandardize them
  scale_attribs_scale_df <- lapply(1:length(scale_attribs_scale), function(i)
    tibble(feat = names(scale_attribs_scale[[i]]),
           value = scale_attribs_scale[[i]],
           LM = i)) %>% 
    do.call("rbind", .) %>% 
    mutate(LM = LM - 1) %>% 
    rename(scale = value)
  
  model_weights <- model_weights %>% 
    left_join(scale_attribs_scale_df, by = join_by(feat, LM)) %>% 
    mutate(coef_unscaled = coef/scale)
  
  var_imp <- var_imp %>% 
    add_row(variable = model_weights$feat,
            value = model_weights$coef,
            train_set = str_replace(f, path_data_complete, ""),
            model = model,
            LM = model_weights$LM,
            var_imp_type = "scaled") %>% 
    add_row(variable = model_weights$feat,
            value = model_weights$coef_unscaled,
            train_set = str_replace(f, path_data_complete, ""),
            model = model,
            LM = model_weights$LM,
            var_imp_type = "unscaled")
  
  # hyperparams are always 0 (not tuned)
  hyperparams <- hyperparams %>%
    add_row(hyperparameter = c("Lam1"),
            value = "0",
            train_set = str_replace(f, path_data_complete, ""),
            model = model,
            LM = 0:max_LM) %>% 
    add_row(hyperparameter = c("Lam2"),
            value = "0",
            train_set = str_replace(f, path_data_complete, ""),
            model = model,
            LM = 0:max_LM) 
  
  # save timings 
  timings <- timings %>% 
    add_row(type = c("tuning", "build final model", "predict"),
            value = c(time_tuning, time_final_model, time_predict),
            train_set = str_replace(f, path_data_complete, ""),
            model = model,
            LM = NA_real_) 
  
  message(sprintf("DONE in %s minutes.", 
                  difftime(Sys.time(), start_time, units = "mins") %>% as.numeric()))
}

# save predictions, variable importance, ...
save(predictions, file = preds_file)
save(var_imp, file = var_imp_file)
save(hyperparams, file = hyperparams_file)
save(timings, file = timings_file)

