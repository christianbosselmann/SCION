# SCN
# Functional variant prediction for voltage-gated sodium channels
# helper script to produce a look-up table for the Brunklaus decision rule
# brunklaus.R

# pkg
library(librarian)
librarian::shelf(tidyverse,
                 readxl, 
                 caret)

# set seed
set.seed(seed)

# get raw data for aa substitution and cid for comparison
dat_raw <- read_csv("data/clean_tbl.csv")
cid <- read_csv("features/cid.csv")[,1:10]

# merge in cid, maintain row order
cid <- cid %>%
  pivot_longer(!cid, names_to = "gene", values_to = "pos") %>%
  group_by(gene) %>%
  arrange(pos) %>%
  distinct(pos, .keep_all = TRUE)

vlookup <- merge(dat_raw, cid) %>%
  arrange(id)

vlookup <- vlookup %>%
  select(-pheno)

# export as csv
# write_csv(dat_raw, "scn_viewer/vlookup.csv")

# set up cv ids
cv <- caret::createFolds(dat_raw$y, k = k)

# initialize loop objects
vlookup_train <- data.frame() 
class_metrics <- metric_set(accuracy, kap, mcc, f_meas)
metrics <- list()

for (i in 1:length(cv)){
  test <- cv[[i]]
  train <- unlist(cv[-i])
  
  # store training subset of Brunklaus lookup table
  vlookup_train <- vlookup[train,]
  
  # count cid per position
  tmp_cid <- vlookup_train %>%
    select(cid, y) 
  
  vlookup_train <- cbind(tmp_cid, model.matrix(~ y - 1, vlookup_train)) %>%
    select(-y)
  
  vlookup_train <- aggregate(vlookup_train[2:3], by = vlookup_train[1], function(v) sum(v))
  
  # implicit encoding, then sum up to get a GOF/LOF score per index position
  # GOF +ve, LOF -ve
  vlookup_train$y_hat <- vlookup_train$yGOF - vlookup_train$yLOF 
  
  # translate into categorical label
  vlookup_train$y_hat[vlookup_train$y_hat > 0] <- "GOF"
  vlookup_train$y_hat[vlookup_train$y_hat < 0] <- "LOF"
  vlookup_train$y_hat[vlookup_train$y_hat == 0] <- "Uncertain"
  
  # optional: drop uncertain labels
  # vlookup_train <- vlookup_train %>%
  #   filter(!y_hat == "Uncertain")
  
  # optional: replace uncertain label with coin toss
  coin <- c("GOF", "LOF")
  vlookup_train$y_hat[vlookup_train$y_hat == "Uncertain"] <- sample(coin, 1)
  
  # use vlookup_train tbl to apply decision rule to test observations
  vlookup_test <- vlookup[test,]
  vlookup_test <- merge(vlookup_train, vlookup_test, by = "cid")
  
  # encode factors
  vlookup_test$y <- factor(vlookup_test$y, levels = c("GOF", "LOF"))
  vlookup_test$y_hat <- factor(vlookup_test$y_hat, levels = c("GOF", "LOF"))
  
  # get metrics
  metrics[[i]] <- vlookup_test %>%
    class_metrics(truth = y, estimate = y_hat, na.rm = FALSE)
}

# create summary and export
report <- rbindlist(metrics, idcol = "fold") %>%
  group_by(.metric) %>%
  summarise(mean = mean(.estimate, na.rm = TRUE), sd = sd(.estimate, na.rm = TRUE))

write_csv(report, paste0('out/report_', Sys.time(), '.csv'))
