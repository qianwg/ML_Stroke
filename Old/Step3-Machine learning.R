rm(list=ls())
library(tidyverse)
library(arrow)
source('functions.R')
# Models
# SVM kernel C gamma
# KNN n_neighbors weights algorithm
# Logistic regression solver penalty C 
# Random forest n_estimators criterion max_depth min_samples_leaf min_samples_split  max_features
# LGBM n_estimators max_depth subsample colsample_bytree learning_rate num_leaves
# XGBoost n_estimators max_depth min_child_weight subsample eta 
# ANN Learning rate Number of layers Layer size Batch size Epochs Dropout optimizer


# Load Data----
DataStep1 = read_feather('D:/UKBiobank/ML_Stroke/Data/memb_map_df_all.feather')
colnames(DataStep1)


# Transform Variable Type----
covarites = colnames(DataStep1)[c(55,56,57,58,64,65,66,
                                  80:120,123,124:153,28:54,
                                  14:27,
                                  28:54,2:9, 156, 159,162,165,168,
                                  171,174,177,180,183,186,189,192,195)]
DataStep2 = DataStep1 %>% select(covarites)
colnames(DataStep2)
DataStep3 = DataStep2 %>%
  mutate_at(c(2,3,5,6,7,24:34,
              35:75,78), as.factor) %>% 
  select(1:76,80:106)

# Select P < 0.05 and missing > 50% in Table 1----
Table1 = read.csv('./tables/V2-Table1.csv')
Table1 %>% filter(p>0.05) %>% 
  select(X)

DataStep4 = DataStep3 %>% 
  select(-c(`Famli history - Bowel cancer`,`Famli history - Breast cancer`, `Famli history - Lung cancer`,
            `Peripheral vasodilators`, `Vasoprotectives`, `Liver disease`, `Renal disease`,
            `Glucocorticoids`, `Dementia`, `Venous thrombosis`:`Rectal cancer`, 
            `Parkinsons disease`, `Cataracts`, `Glaucoma`,
            `Red blood cell (erythrocyte) distribution width`, `Mean platelet (thrombocyte) volume`,
            `C-reactive protein`, `Glycated haemoglobin (HbA1c)`, `Lipoprotein A`))

# Missing values----
colnames(DataStep4) <- gsub(" ", "_", colnames(DataStep4))

# colnames(DataStep4)[c(6,11,12,16,17,18,19:24)] = c('Body_mass_index', 'Systolic_blood_pressure',
#                                     'Diastolic_blood_pressure','FVC',
#                                     'FEV1', 'FEV1_FVC_ratio', 'Famliy_history_Dementia',
#                                     'Famliy_history_Chronic_bronchitis_and_emphysema',
#                                     'Famliy_history_Heart_disease',
#                                     'Famliy_history_High_blood_pressure',
#                                     'Famliy_history_Severe_depression',
#                                     'Famliy_history_Stroke')

colnames(DataStep4)[c(8,13,14,15,21:31,75)] = c('Body_mass_index','Systolic_blood_pressure',
                                        'Diastolic_blood_pressure', 'Pulse_rate', 'FVC','FEV1', 'FEV1_FVC_ratio',
                                        'Famli_history_Dementia', 'Famli_history_Chronic_bronchitis_and_emphysema',
                                        'Famli_history_Diabetes', 'Famli_history_Heart_disease',
                                        'Famli_history_High_blood_pressure', 'Famli_history_Parkinsons_disease',
                                        'Famli_history_Severe_depression', 'Famli_history_Stroke',
                                        'White_blood_cell')


missing_perc_cols <- colMeans(is.na(DataStep4))  # 计算每列的缺失值百分比
df_filtered_cols <- DataStep4[, missing_perc_cols <= 0.25]  # 保留缺失值少于 25% 的列

missing_perc_rows <- rowMeans(is.na(df_filtered_cols))  # 计算每行的缺失值百分比
df_filtered_rows <- df_filtered_cols[missing_perc_rows <= 0.30, ]  # 保留缺失值少于 30% 的行

DataStep5 = df_filtered_rows
col_missing = DataStep5 %>% 
  select_if(~any(is.na(.))) %>% 
  colnames()
set.seed(1234)
imp <- mice(DataStep5 %>% select(col_missing))
Data_missing <- complete(imp, action = 3) %>% tibble()

DataStep6 = DataStep5
DataStep6[col_missing] = Data_missing

# saveRDS(DataStep6, './Data/DataStep6.rds')

DataStep6 = readRDS('./Data/DataStep6.rds')


# information gain
library(FSelector)
ig <- information.gain(Cluster ~ ., data = DataStep6)
ranked_features <- ig[order(-ig$attr_importance), , drop = FALSE]
print(ranked_features)
ranked_features %>% write.csv('tables/importance.csv')
selected_features <- rownames(ranked_features)[1:10]

DataStep7 = DataStep6[c('Cluster', selected_features)]
# saveRDS(DataStep7, './Data/V2-DataStep7.rds')
# DataStep7 %>% write_feather('./Data/V2-DataStep7.feather')


DataStep7 = readRDS('./Data/V2-DataStep7.rds')


# Tidymodels-recipes----
library(tidymodels)
set.seed(222)
data_split <- initial_split(DataStep7, prop = 3/4, strata = 'Cluster')

train_data <- training(data_split)
test_data  <- testing(data_split)

Pheno_rec <- 
  recipe(Cluster ~ ., data = train_data) %>% 
  update_role(eid, new_role = "ID") %>% 
  step_normalize(all_numeric()) %>%
  step_dummy(all_nominal_predictors(), one_hot = TRUE) %>%
  step_zv()

Pheno_prep <- prep(Pheno_rec)
juiced <- juice(Pheno_prep)


set.seed(1234)
Pheno_folds <- vfold_cv(train_data, v = 10)

library(future)
plan(multisession)

all_cores <- parallel::detectCores(logical = FALSE)

library(doParallel)
cl <- makePSOCKcluster(all_cores)
registerDoParallel(cl)


# Random forest ----
rf_spec <- rand_forest(
  mtry = tune(),
  trees = tune(),
  min_n = tune()
) %>% 
  set_mode("classification") %>% 
  set_engine("ranger")

rf_wf <- workflow() %>% 
  add_recipe(Pheno_rec) %>% 
  add_model(rf_spec)


rf_params <- extract_parameter_set_dials(rf_spec) %>% 
  update(mtry = mtry(range = c(1, 10)))

rf_grid = grid_space_filling(
  rf_params,
  size = 100
)

set.seed(456)
rf_tune_res <- tune_grid(
  rf_wf,
  resamples = Pheno_folds,
  grid = rf_grid
)

rf_tune_res %>%
  collect_metrics() %>%
  filter(.metric == "roc_auc") %>%
  select(mean, mtry:min_n) %>%
  pivot_longer(mtry:min_n,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(alpha = 0.8, show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "AUC")

saveRDS(rf_tune_res, './Data/RF_regular_res.rds')

best_auc <- select_best(rf_tune_res,metric = 'roc_auc')

final_rf <- finalize_model(
  rf_spec,
  best_auc
)


final_wf <- workflow() %>%
  add_recipe(Pheno_rec) %>%
  add_model(final_rf)

final_res <- final_wf %>%
  last_fit(data_split)

final_res %>%
  collect_metrics()

final_res %>%
  collect_predictions() %>%
  roc_curve(Cluster, .pred_1) %>% 
  autoplot()

# Xgboost ----
xgb_spec <- boost_tree(
  trees = tune(),
  tree_depth = tune(), min_n = tune(),
  loss_reduction = tune(),                     ## first three: model complexity
  sample_size = tune(), mtry = tune(),         ## randomness
  learn_rate = tune()                          ## step size
) %>% 
  set_mode("classification") %>% 
  set_engine("xgboost")


xbg_wf <- workflow() %>% 
  add_recipe(Pheno_rec) %>% 
  add_model(xgb_spec)
  
xgb_grid <- grid_space_filling(
  tree_depth(),
  min_n(),
  loss_reduction(),
  sample_size = sample_prop(),
  finalize(mtry(), train_data),
  learn_rate(),
  trees(),
  size = 500
)

set.seed(456)
xbg_tune_res <- tune_grid(
  xbg_wf,
  resamples = Pheno_folds,
  grid = xgb_grid
)

saveRDS(xbg_tune_res, './Data/xbg_tune_res.rds')


xbg_tune_res %>%
  collect_metrics() %>%
  filter(.metric == "roc_auc") %>%
  filter(mean>0.84) %>% 
  select(mean, mtry:sample_size) %>%
  pivot_longer(mtry:sample_size,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(alpha = 0.8, show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "AUC")

show_best(xbg_tune_res,metric = "roc_auc")
best_auc <- select_best(xbg_tune_res,metric= "roc_auc")
best_auc


final_xgb <- finalize_model(
  xgb_spec,
  best_auc
)


final_wf <- workflow() %>%
  add_recipe(Pheno_rec) %>%
  add_model(final_xgb)

final_res <- final_wf %>%
  last_fit(data_split)

final_res %>%
  collect_metrics()

final_res %>%
  collect_predictions() %>%
  roc_curve(Cluster, .pred_1) %>% 
  autoplot()


# KNN----
KNN_spec = nearest_neighbor(
  neighbors = tune(),
  weight_func = tune()
) %>% 
  set_engine("kknn") %>% 
  set_mode("classification")

KNN_wf <- workflow() %>% 
  add_recipe(Pheno_rec) %>% 
  add_model(KNN_spec)

# KNN_grid <- expand_grid(
#   neighbors = seq(10,100,10),
#   weight_func = c("rectangular", "triangular", "epanechnikov", "biweight", 
#               "triweight", "cos", "inv", "gaussian", "rank", "optimal")
# )

KNN_grid <- expand_grid(
  neighbors = seq(250,350,10),
  weight_func = c("rectangular", "triangular", "epanechnikov",
                  "cos", "inv", "gaussian", "rank", "optimal"),
)

KNN_params <- extract_parameter_set_dials(KNN_spec) %>%
  update(
    neighbors = neighbors(range = c(700, 1000)),
    weight_func = weight_func(values = c("triangular", "epanechnikov",
                                         "cos"))
  )

KNN_grid = grid_space_filling(KNN_params,
                              size = 100)

set.seed(456)
KNN_tune_res <- tune_grid(
  KNN_wf,
  resamples = Pheno_folds,
  grid = KNN_grid
)

KNN_tune_res %>%
  collect_metrics() %>%
  filter(.metric == "roc_auc") %>%
  mutate(weight_func = factor(weight_func)) %>%
  ggplot(aes(neighbors, mean, color = weight_func)) +
  geom_line(alpha = 0.5, size = 1.5) +
  geom_point() +
  labs(y = "AUC")

saveRDS(KNN_tune_res, './Data/KNN_tune_res.rds')
show_best(KNN_tune_res,metric = "roc_auc")
best_auc <- select_best(KNN_tune_res,metric= "roc_auc")
best_auc

final_KNN <- finalize_model(
  KNN_spec,
  best_auc
)


final_wf <- workflow() %>%
  add_recipe(Pheno_rec) %>%
  add_model(final_KNN)

final_res <- final_wf %>%
  last_fit(data_split)

final_res %>%
  collect_metrics()

final_res %>%
  collect_predictions() %>%
  roc_curve(Cluster, .pred_1) %>% 
  autoplot()

# SVM linear----
SVM_linear_spec = svm_linear(cost = tune()) %>% 
  set_mode("classification") %>% 
  set_engine("kernlab")

SVM_linear_wf <- workflow() %>% 
  add_recipe(Pheno_rec) %>% 
  add_model(SVM_linear_spec)

SVM_linear_grid <- expand_grid(
  cost = 10^seq(log10(0.0001), log10(10000))
)

set.seed(456)
SVM_linear_tune_res <- tune_grid(
  SVM_linear_wf,
  resamples = Pheno_folds,
  grid = SVM_linear_grid
)

SVM_linear_tune_res %>%
  collect_metrics() %>%
  filter(.metric == "roc_auc") %>%
  ggplot(aes(cost, mean)) +
  geom_line(alpha = 0.5, size = 1.5) +
  geom_point() +
  labs(y = "AUC")

saveRDS(SVM_linear_tune_res, './Data/SVM_linear_tune_res.rds')
show_best(SVM_linear_tune_res,metric = "roc_auc")
best_auc <- select_best(SVM_linear_tune_res,metric= "roc_auc")
best_auc

final_SVM_linear <- finalize_model(
  SVM_linear_spec,
  best_auc
)


final_wf <- workflow() %>%
  add_recipe(Pheno_rec) %>%
  add_model(final_SVM_linear)

final_res <- final_wf %>%
  last_fit(data_split)

final_res %>%
  collect_metrics()

final_res %>%
  collect_predictions() %>%
  roc_curve(Cluster, .pred_1) %>% 
  autoplot()

# SVM ploy----





# SVM rbf----
SVM_rbf_spec = svm_rbf(
  cost = tune(),
  rbf_sigma = tune()) %>% 
  set_mode("classification") %>% 
  set_engine("kernlab")

SVM_rbf_wf <- workflow() %>% 
  add_recipe(Pheno_rec) %>% 
  add_model(SVM_rbf_spec)


SVM_rbf_params <- extract_parameter_set_dials(SVM_rbf_spec)

SVM_rbf_grid <- grid_space_filling(
  SVM_rbf_params,
  size = 30
)

set.seed(456)
SVM_rbf_tune_res <- tune_grid(
  SVM_rbf_wf,
  resamples = Pheno_folds,
  grid = SVM_rbf_grid
)

saveRDS(SVM_rbf_tune_res, './Data/SVM_rbf_tune_res.rds')

SVM_rbf_tune_res %>%
  collect_metrics() %>%
  filter(mean>0.8) %>%
  select(mean, cost:rbf_sigma) %>%
  pivot_longer(cost:rbf_sigma,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_line(alpha = 0.8, show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "AUC")


show_best(SVM_rbf_tune_res,metric = "roc_auc")
best_auc <- select_best(SVM_rbf_tune_res,metric= "roc_auc")
best_auc

final_SVM_rbf <- finalize_model(
  SVM_rbf_spec,
  best_auc
)

final_wf <- workflow() %>%
  add_recipe(Pheno_rec) %>%
  add_model(final_SVM_rbf)

final_res <- final_wf %>%
  last_fit(data_split)

final_res %>%
  collect_metrics()

final_res %>%
  collect_predictions() %>%
  roc_curve(Cluster, .pred_1) %>% 
  autoplot()


# logistic----
logistic_spec = logistic_reg(
  penalty = tune(),
  mixture = tune()
) %>%
  set_engine("glmnet") %>%
  set_mode("classification") 

logistic_wf <- workflow() %>% 
  add_recipe(Pheno_rec) %>% 
  add_model(logistic_spec)

logistic_params <- extract_parameter_set_dials(logistic_spec)

logistic_grid <- grid_space_filling(
  logistic_params,
  size = 100
)


set.seed(456)
logistic_tune_res <- tune_grid(
  logistic_wf,
  resamples = Pheno_folds,
  grid = logistic_grid
)

saveRDS(logistic_tune_res, './Data/logistic_tune_res.rds')

logistic_tune_res %>%
  collect_metrics() %>%
  # mutate(penalty = as.factor(penalty)) %>% 
  filter(.metric == "roc_auc") %>%
  ggplot(aes(penalty, mean)) +
  geom_line(alpha = 0.5, size = 1.5) +
  geom_point() +
  labs(y = "AUC")

show_best(logistic_tune_res,metric = "roc_auc")
best_auc <- select_best(logistic_tune_res,metric= "roc_auc")
best_auc

final_logistic_rbf <- finalize_model(
  logistic_spec,
  best_auc
)

final_wf <- workflow() %>%
  add_recipe(Pheno_rec) %>%
  add_model(final_logistic_rbf)

final_res <- final_wf %>%
  last_fit(data_split)

final_res %>%
  collect_metrics()

final_res %>%
  collect_predictions() %>%
  roc_curve(Cluster, .pred_1) %>% 
  autoplot()


# LGBM----
LGBM_spec = boost_tree(
  mtry = tune(), 
  trees = tune(), 
  tree_depth = tune(), 
  learn_rate = tune(), 
  min_n = tune(), 
  loss_reduction = tune()
) %>%
  set_engine("lightgbm") %>%
  set_mode("classification")

library(bonsai)
LGBM_wf <- workflow() %>% 
  add_recipe(Pheno_rec) %>% 
  add_model(LGBM_spec)

LGBM_params <- extract_parameter_set_dials(LGBM_spec) %>% 
  update(mtry = mtry(range = c(1, 10)))

LGBM_grid <- grid_space_filling(
  LGBM_params,
  size = 300
)

set.seed(456)
LGBM_tune_res <- tune_grid(
  LGBM_wf,
  resamples = Pheno_folds,
  grid = LGBM_grid
)

saveRDS(LGBM_tune_res , './Data/LGBM_tune_res.rds')
LGBM_tune_res = readRDS('./Data/LGBM_tune_res.rds')

LGBM_tune_res %>%
  collect_metrics() %>%
  filter(.metric == "roc_auc") %>%
  select(mean, mtry:loss_reduction) %>%
  pivot_longer(mtry:loss_reduction,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(alpha = 0.8, show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "AUC")

show_best(LGBM_tune_res,metric = "roc_auc")
best_auc <- select_best(LGBM_tune_res,metric= "roc_auc")
best_auc

final_LGBM_rbf <- finalize_model(
  LGBM_spec,
  best_auc
)

final_wf <- workflow() %>%
  add_recipe(Pheno_rec) %>%
  add_model(final_LGBM_rbf)

final_res <- final_wf %>%
  last_fit(data_split)

final_res %>%
  collect_metrics()

final_res %>%
  collect_predictions() %>%
  roc_curve(Cluster, .pred_1) %>% 
  autoplot()


# ANN----
nnet_spec <-
  mlp(hidden_units = tune(), penalty = tune(), epochs = tune()) %>%
  set_mode("classification") %>%
  set_engine("brulee")



# GAM----


# decision tree----


