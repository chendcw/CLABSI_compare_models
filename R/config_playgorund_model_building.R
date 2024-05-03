# libraries
library(glmnet)
library(doParallel)

registerDoParallel(detectCores())

# common config
# -------------

N_iter <- 100

frac_test <- 1/3

# models_path
# -----------

models_path <- "models/"
var_sel_path <- "playground/variable_selection/selection_by_reviewers/"

# BASELINE config
# ---------------

data_path_play_base_miss <- "data_for_models/2012_2013/missing/BASE/"
data_path_play_base_complete_MF <- "data_for_models/2012_2013/imputed/BASE/MF/"

cat_cols_base <- c("MS_medical_specialty", 
                   "GEN_LM_day_categ", 
                   "GEN_LM_month_categ", 
                   "GEN_LM_season_categ")

var_imp_dir_base <- "variable_importance/2012_2013/BASE/"
preds_dir_base <- "predictions/2012_2013/BASE/"



vars_selected_DOMK_LIT <- c("CAT_catheter_type_binary_all_CVC", 
                            "CAT_catheter_type_binary_all_Port_a_cath",
                            "CAT_catheter_type_binary_all_Tunneled_CVC", 
                            "CAT_catheter_type_binary_all_PICC", 
                            "CAT_catheter_location_binary_all_Collarbone",
                            "CAT_catheter_location_binary_all_Neck",
                            # "MED_TPN",
                            "MED_7d_TPN", # name change
                            # "MED_L2_J01_ANTIBACTERIALS_FOR_SYSTEMIC_USE",
                            "MED_L2_7d_J01_ANTIBACTERIALS_FOR_SYSTEMIC_USE", # name change
                            # "MED_L2_L01_ANTINEOPLASTIC_AGENTS",
                            "MED_L2_7d_L01_ANTINEOPLASTIC_AGENTS", # name change
                            "CLABSI_history",
                            "COM_PATH_tumor_before_LM",
                            "CARE_VS_temperature_max",
                            "CARE_VS_systolic_BP_last",
                            "LAB_WBC_count_last",
                            # "LAB_WBC_Neutrophils_last", # decided to exclude
                            # "CARE_VS_MV",
                            "COM_lymphoma_before_LM",
                            "COM_PATH_transplant_before_LM",
                            #"MB_other_infection_than_BSI", # decided to exclude
                            "MB_other_infection_than_BSI_during_window",
                            "LAB_CRP_last",
                            "ADM_admission_source_binary_all_Home",
                            "CARE_VS_MV",
                            "MS_is_ICU_unit")


# DYNAMIC config
# --------------
data_path_play_dyn_miss <- "data_for_models/2012_2013/missing/DYN/"
data_path_play_dyn_complete_MF <- "data_for_models/2012_2013/imputed/DYN/MF/"


