# SCN
# Functional variant prediction for voltage-gated sodium channels
# helper script to provide a baseline comparison to the Heyne GBM
# code adapted from github.com/heyhen/fuNCion
# helper_heyne.R

# pkg
library(librarian)
librarian::shelf(tidyverse,
                 caret)

# set seed
set.seed(42)

# get raw data for aa substitution and cid for comparison
data <- read_csv("data/dat_prep.csv")

# 
fitControl <- caret::trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 10, 
  classProbs = TRUE
  )

inTraining <- createDataPartition(as.factor(data$y), p = .9, list = FALSE) 
training <- data[ inTraining,] 
testing <- data[ -inTraining,] 

# updated upsampling procedure from Heyne code
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

modelperformance <- function(out) {
  res <- c(multiClassSummary(out, lev = c("GOF", "LOF")),
           mcc(out, obs, pred)$.estimate, # updated from Heyne code due to outdated dependency
           round(twoClassSummary(out, lev = c("GOF", "LOF")), digits = 2))
  names(res)[15] <- "MCC"
  return(res[c("Balanced_Accuracy", "Sens", "Spec","AUC","Precision","Recall","F1", "prAUC","Kappa", "MCC")])
}

# single 0.1 holdout as in original Heyne method, then create report
report_metrics <- as.data.frame(modelperformance(out)) %>%
  rename(value = 'modelperformance(out)')
report_metrics$metric <- rownames(report_metrics) 
rownames(report_metrics) <- NULL
report_metrics <- report_metrics %>%
  select(metric, everything()) 


# export output
write_csv(report_metrics, paste0('out/report_metrics_', Sys.time(), '.csv'))
write_csv(out, paste0('out/report_preds_', Sys.time(), '.csv'))
write_csv(model1$bestTune, paste0('out/report_params_', Sys.time(), '.csv'))
