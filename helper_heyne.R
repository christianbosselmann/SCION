# SCN
# Functional variant prediction for voltage-gated sodium channels
# helper script to provide a baseline comparison to the Heyne GBM
# code adapted from github.com/heyhen/fuNCion
# helper_heyne.R

# pkg
library(librarian)
librarian::shelf(tidyverse,
                 tidymodels,
                 caret,
                 data.table)

# set seed
set.seed(seed)

# get data
data <- read_csv("data/dat_prep.csv") 
data <- data %>% select(-gene) # task, not a feature

# set splits
cv <- caret::createFolds(data$y, k = k)

# set fit control
fitControl <- caret::trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 10, 
  classProbs = TRUE
)

# k-fold cv on same splits (seed)
out <- list()
for (i in 1:length(cv)){
  ind_test <- cv[[i]]
  ind_train <- unlist(cv[-i], use.names = F)
  
  testing <- data[ind_test,] 
  training <- data[ind_train,] 

  # note: no additional feature engineering as in original Heyne code,
  # due to the different data set
  model <- caret::train(y ~ ., data = training,
                         method = "gbm",
                         trControl = fitControl,
                         verbose = FALSE)
  
  test_data <- testing$y
  
  out[[i]] <- data.frame(obs = as.factor(test_data),
                    GOF = predict(model, newdata = testing, type = "prob")[,"GOF"],
                    LOF = predict(model, newdata = testing, type = "prob")[,"LOF"],
                    pred = predict(model, newdata = testing),
                    fold = i)
}
out <- rbindlist(out)

# create report
class_metrics <- metric_set(accuracy, kap, mcc, sens, spec, f_meas, roc_auc, pr_auc)

report_metrics <- out %>%
  group_by(fold) %>%
  class_metrics(truth = obs, GOF, estimate = pred) %>% 
  group_by(.metric) %>%
  summarise(mean = mean(.estimate, na.rm = TRUE), sd = sd(.estimate, na.rm = TRUE))

# export output
write_csv(report_metrics, paste0('out/report_metrics_', Sys.time(), '.csv'))
write_csv(out, paste0('out/report_preds_', Sys.time(), '.csv'))
write_csv(model$bestTune, paste0('out/report_params_', Sys.time(), '.csv'))
