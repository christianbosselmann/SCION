# SCN
# Functional variant prediction for voltage-gated sodium channels
# helper script to provide a baseline comparison to the Heyne GBM
# code adapted from github.com/heyhen/fuNCion
# helper_heyne.R

# pkg
library(librarian)
librarian::shelf(tidyverse,
                 tidymodels,
                 caret)

# set seed
set.seed(seed)

# get data
data <- read_csv("data/dat_prep.csv") 
data <- data %>% select(-gene) # task, not a feature

# set fit control
fitControl <- caret::trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 10, 
  classProbs = TRUE
  )

inTraining <- createDataPartition(as.factor(data$y), p = .9, list = FALSE) 
training <- data[ inTraining,] 
testing <- data[ -inTraining,] 

# upsampling, as in Heyne
training <- upSample(x = training[, -ncol(training)],
                     y = as.factor(training$y)) 

training <- training %>% rename(y = Class)

# note: no additional feature engineering as in original Heyne code,
# due to the different data set
model1 <- caret::train(y ~ ., data = training,
                       method = "gbm",
                       trControl = fitControl,
                       verbose = FALSE)

test_data <- testing$y

out <- data.frame(obs = as.factor(test_data),
                  GOF = predict(model1, newdata = testing, type = "prob")[,"GOF"],
                  LOF = predict(model1, newdata = testing, type = "prob")[,"LOF"],
                  pred = predict(model1, newdata = testing)
)

# single 0.1 holdout as in original Heyne method, then create report
class_metrics <- metric_set(accuracy, kap, mcc, sens, spec, f_meas, roc_auc, pr_auc)

report_metrics <- out %>%
  class_metrics(truth = obs, GOF, estimate = pred) 

# export output
write_csv(report_metrics, paste0('out/report_metrics_', Sys.time(), '.csv'))
write_csv(out, paste0('out/report_preds_', Sys.time(), '.csv'))
write_csv(model1$bestTune, paste0('out/report_params_', Sys.time(), '.csv'))
