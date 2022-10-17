if (!require("pacman")) {
  install.packages("pacman")
}

pacman::p_load(
  tidyverse,
  zeallot,
  tidymodels,
  xgboost,
  mlrMBO,
  SHAPforxgboost
)

# global ------------------------------------------------------------------

seed <- 1234
set.seed(seed)

tree_method <- "hist"
temp_path <- "./data/temp/"

save_model <- T

n_iter <- 60

# data --------------------------------------------------------------------

data_raw <- read_csv("./data/events_6730330.csv")

# data quality check

data_raw %>%
  sapply(function(x) sum(is.na(x))/length(x))

# Preprocess --------------------------------------------------------------

data_process <- data_raw %>%
  select(
    flood_peak,
    rain_duration_time:flood_type
  ) %>%
  mutate(
    id = 1:n()
  )

data_recipe <- recipe(x = data_process) %>%
  update_role(rain_duration_time:flood_type, new_role = "predictor") %>%
  update_role(flood_peak, new_role = "outcome") %>%
  update_role(id, new_role = "id variable") %>%
  step_dummy(all_nominal_predictors(), one_hot = T) %>%
  prep()

data_process <- bake(data_recipe, new_data = data_process) 

predictor_names <- data_process %>%
  names() %>%
  setdiff(c("flood_peak", "id"))

# resampling -------------------------------------------------------------

cv_splits <- rsample::nested_cv(
  data_process, 
  inside = vfold_cv(v = 5, strata = flood_peak),
  outside = vfold_cv(v = 5, strata = flood_peak)
)


# Modeling ----------------------------------------------------------------

# outer CV iteration
preds <- vector("list", 5)
SHAPs <- vector("list", 5)
SHAP_longs <- vector("list", 5)

for (i in 1:5) {
  
  outer_cv_iter <- i
  
  outer_split <- cv_splits$splits[[outer_cv_iter]]
  inner_split <- cv_splits$inner_resamples[[outer_cv_iter]]
  
  # prepare data sets
  training_outer <- analysis(outer_split)
  training_outer_id <- training_outer %>%
    pull(id) %>%
    unlist()
  
  test <- assessment(outer_split)
  test_id <- test %>%
    pull(id) %>%
    unlist()
  
  val <- lapply(inner_split$splits, assessment)
  val <- val %>% lapply(function(x) x %>% pull(id) %>% unlist())
  val <- lapply(val, function(x)
    which(training_outer_id %in% x)) # get relative location
  
  # prepare datasets
  dtrain <- xgb.DMatrix(data = data.matrix(training_outer[predictor_names]),
                        label = training_outer[outcome_names(data_recipe)] %>% unlist() %>% unname())
  
  dtest <- xgb.DMatrix(data = data.matrix(test[predictor_names]),
                       label = test[outcome_names(data_recipe)] %>% unlist() %>% unname())
  
  dall <- xgb.DMatrix(data = data.matrix(data_process[predictor_names]),
                      label = data_process[outcome_names(data_recipe)] %>% unlist() %>% unname())
  
  # scoring function
  scoringFunction <- function(x){
    eta <- x["eta"] %>% unlist()
    max_depth <- x["max_depth"] %>% unlist()
    min_child_weight <- x["min_child_weight"] %>% unlist()
    subsample <- x["subsample"] %>% unlist()
    colsample_bytree <- x["colsample_bytree"] %>% unlist()
    gamma <- x["gamma"] %>% unlist()
    
    Pars <- list(
      booster = "gbtree",
      eta = eta,
      max_depth = max_depth,
      min_child_weight = min_child_weight,
      colsample_bytree = colsample_bytree,
      subsample = subsample,
      gamma = gamma
    )
    
    xgbcv <- xgb.cv(
      objective = "reg:squarederror",
      data = dtrain,
      folds = val,
      tree_method = tree_method,
      max_bin = 256,
      nround = 5000,
      early_stopping_rounds = 20,
      verbose = 0,
      params = Pars,
    )
    
    best_iteration <- xgbcv$best_iteration
    watchlist <- NULL
    xgbFit <- xgb.train(
      data = dtrain,
      objective = "reg:squarederror",
      tree_method = tree_method,
      max_bin = 256,
      nround = best_iteration,
      verbose = 0,
      params = Pars
    )
    
    # save model
    if (save_model){
      
      model_id <- (list.files(temp_path, pattern = "\\.model$") %>% length()) + 1
      file_path <- paste0(temp_path, "model_",model_id,".model")
      xgb.save(model = xgbFit, fname = file_path)
    }
    
    # output
    min(xgbcv$evaluation_log$test_rmse_mean)
  }
  
  obj_fun <- makeSingleObjectiveFunction(
    fn = scoringFunction,
    par.set = makeParamSet(
      makeNumericParam("eta",                    lower = 0.005, upper = 0.1),
      makeIntegerParam("max_depth",              lower= 2,      upper = 10),
      makeIntegerParam("min_child_weight",       lower= 1,    upper = 10),
      makeNumericParam("subsample",              lower = 0.20,  upper = 1),
      makeNumericParam("colsample_bytree",       lower = 0.20,  upper = 1),
      makeNumericParam("gamma",                  lower = 0,     upper = 10)
    ),
    has.simple.signature = FALSE,
    minimize = TRUE
  )
  
  # run optimization
  
  unlink(paste0(temp_path, "/*")) # clean model path
  
  des = generateDesign(
    n = 4 * getNumberOfParameters(obj_fun),
    par.set = getParamSet(obj_fun),
    fun = lhs::randomLHS
  )
  
  des$y = apply(des, 1, obj_fun)
  
  control <- makeMBOControl() %>%
    setMBOControlTermination(., iters = n_iter - 4 * getNumberOfParameters(obj_fun))
  
  run <- mbo(
    fun = obj_fun,
    design = des,
    control = control,
    show.info = TRUE
  )
  
  # Post-processing ---------------------------------------------------------
  
  xgbFit <- xgboost::xgb.load(paste0(temp_path, "model_",run$best.ind,".model"))
  
  preds[[i]] <- tibble(
    id = test$id,
    pred = predict(xgbFit, dtest),
    ob = getinfo(dtest, "label"),
    iter = i
  )
  
  dataX <- data.matrix(data_process[predictor_names])
  shap_values <- shap.values(xgb_model = xgbFit, X_train = dataX)
  SHAPs[[i]] <- shap_values
  
  shap_long <- shap.prep(xgb_model = xgbFit, X_train = dataX)
  SHAP_longs[[i]] <- shap.prep(shap_contrib = shap_values$shap_score, X_train = dataX)
}



# Save data ---------------------------------------------------------------

save(preds, SHAP_longs, SHAPs, file = "small_AB.Rda")



# Post-process ------------------------------------------------------------


data_plot <- preds %>%
  bind_rows() %>%
  mutate(exp = paste0("experiment #", iter))

data_plot %>%
  ggplot(aes(pred, ob, color = exp, shape = exp)) +
  geom_point() +
  labs(x = "Predicted flood peak",
       y = "Observed flood peak") +
  theme_bw()

hydroGOF::gof(data_plot$pred, data_plot$ob)






load("small_AB.Rda")

outs <- vector("list", 5)

for (i in 1:5){
  outs[[i]] <- tibble(feature = SHAPs[[i]]$mean_shap_score %>% names(),
                      importance =  SHAPs[[i]]$mean_shap_score,
                      exp = paste0("experiment ", i))
}

data_plot <- outs %>%
  bind_rows()

feature_level <- data_plot %>%
  group_by(feature) %>%
  summarise(mean_importance = mean(importance)) %>%
  arrange(mean_importance) %>%
  pull(feature)

data_plot <- data_plot %>%
  mutate(feature = factor(feature, levels = feature_level))

ggplot(data_plot, aes(feature, importance, color = exp, shape = exp)) +
  geom_point()+
  geom_line(aes(group = exp, linetype = exp)) +
  coord_flip()+
  theme_bw()



data_plot <- outs %>%
  bind_rows()

data_plot2 <- data_plot %>%
  group_by(exp) %>%
  mutate(importance_rank = rank(-importance)) %>%
  mutate(feature = factor(feature, levels = feature_level))


ggplot(data_plot2, aes(feature, importance_rank, color = exp, shape = exp)) +
  geom_point() +
  geom_line(aes(group = exp, linetype = exp)) +
  coord_flip()+
  theme_bw()




dall <- xgb.DMatrix(data = data.matrix(data_process[predictor_names]),
                    label = data_process[outcome_names(data_recipe)] %>% unlist() %>% unname())
org <- predict(xgbFit, dall)

data_process2 <- data_process %>%
  mutate(chemistry_Polyamide = 1)

dall <- xgb.DMatrix(data = data.matrix(data_process2[predictor_names]),
                    label = data_process2[outcome_names(data_recipe)] %>% unlist() %>% unname())
mod <- predict(xgbFit, dall)


plot(mod[org!=mod] - org[org!=mod])

tibble(mod = mod,
       org = org,
       id = data_process$id) %>%
  mutate(change = mod - org) %>%
  filter(change != 0) %>%
  ggplot() +
  geom_point(aes(id, change))+
  theme_bw()+
  labs(y = "changes in log10(A/B)")


hydroGOF::gof(data_plot$pred, data_plot$ob)


predictor_names

