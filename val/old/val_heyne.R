# get packages
library(librarian)
librarian::shelf(tidyverse,
                 janitor)

# get data from SCN model output and Heyne
dat_scn <- read.xlsx("app/output.xlsx")
dat_heyne <- read.xlsx("val/heyne_raw.xlsx")

# method:
# Table S4 from Heyne et al. (electrophysiological experiments), filtered for 
# complete cases (NA) and variants not used in model training (in_func_pred = FALSE)
# n = 38

# filter heyne for appropriate variants
dat_heyne <- dat_heyne %>%
  filter(in_func_pred == TRUE) %>%
  filter(pathpred == "pathogenic") %>%
  select(gene, aa1, pos, aa2, obs, heyne_pred, heyne_GOF, heyne_LOF) %>%
  filter(obs == "lof" | obs == "gof") %>%
  filter(complete.cases(.))

dat_heyne$obs <- toupper(dat_heyne$obs)
dat_heyne$heyne_pred <- toupper(dat_heyne$heyne_pred)

# join with our own prediction output
dat_all <- inner_join(dat_heyne, dat_scn)

# fix data types
dat_all <- hablar::retype(dat_all)
dat_all$obs <- as.factor(dat_all$obs)
dat_all$scn_pred <- as.factor(dat_all$scn_pred)
dat_all$heyne_pred <- as.factor(dat_all$heyne_pred)

# get metrics
class_metrics <- metric_set(accuracy, kap, mcc, f_meas, roc_auc, pr_auc)

report_scn <- dat_all %>%
  class_metrics(truth = obs, scn_GOF, estimate = scn_pred)

report_heyne <- dat_all %>%
  class_metrics(truth = obs, heyne_GOF, estimate = heyne_pred)


