# SCN
# Functional variant prediction for voltage-gated sodium channels
# helper script to provide a baseline comparison to the Brunklaus decision rule
# helper_brunklaus.R

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
vlookup_test <- data.frame() 
class_metrics <- metric_set(accuracy, spec, sens, kap, mcc, f_meas)
metrics <- list()
preds <- list()

for (i in 1:length(cv)){
  test <- cv[[i]]
  train <- unlist(cv[-i], use.names = F)
  
  # create subsets
  vlookup_train <- vlookup[train,]
  vlookup_test <- vlookup[test,]
  
  # join in Brunklaus' "mismatch pairs"
  vlookup_test <- left_join(x = vlookup_test, 
                            y = vlookup_train, 
                            by = c("cid", "aa1", "aa2"))
  
  # clean up
  vlookup_test <- vlookup_test %>%
    select(id.x, aa1, gene.x, aa2, pos.x, cid, y.x, y.y) %>%
    rename(id = id.x, gene = gene.x, pos = pos.x, y = y.x, y_hat = y.y)
  
  # all variants not captured by Brunklaus pairs are labeled uncertain
  vlookup_test$y_hat[is.na(vlookup_test$y_hat)] <- "Uncertain"
  
  # encode factors
  vlookup_test$y <- factor(vlookup_test$y, levels = c("GOF", "LOF", "Uncertain"))
  vlookup_test$y_hat <- factor(vlookup_test$y_hat, levels = c("GOF", "LOF", "Uncertain"))

  # get metrics
  metrics[[i]] <- vlookup_test %>%
    class_metrics(truth = y, estimate = y_hat)

  preds[[i]] <- vlookup_test
}

# create summary and export
report_metrics <- rbindlist(metrics, idcol = "fold") %>%
  group_by(.metric) %>%
  summarise(mean = mean(.estimate, na.rm = TRUE), sd = sd(.estimate, na.rm = TRUE))

report_preds <- rbindlist(preds, idcol = "fold")

write_csv(report_metrics, paste0('out/report_metrics_', Sys.time(), '.csv'))
write_csv(report_preds, paste0('out/report_preds_', Sys.time(), '.csv'))
