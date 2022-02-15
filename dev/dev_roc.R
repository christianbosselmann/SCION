# SCN dev
# this script takes "report_final.csv" files, i.e. the raw prediction from nested k-fold cv, and prepares a grid of ROC plots
# just a code snippet, remove before publication

# pkg
library(librarian)
librarian::shelf(ggplot2,
                 ggpubr,
                 grid,
                 pROC)

# data
report_dirac <- read_csv("out/preds dirac.csv")
report_union <- read_csv("out/preds union.csv")
report_mtl <- read_csv("out/preds mtl.csv")
report_mkl <- read_csv("out/preds jaccard semkl.csv")

# ROC
roc_mkl <- report_mkl %>% 
  roc(response = "truth", predictor = "GOF") %>%
  ggroc() +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") +
  theme_bw() +
  coord_fixed() +
  ggtitle("MKL")

roc_dirac <- report_dirac %>% 
  roc(response = "truth", predictor = "GOF") %>%
  ggroc() +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") +
  theme_bw() +
  coord_fixed() +
  ggtitle("Dirac")

roc_union <- report_union %>% 
  roc(response = "truth", predictor = "GOF") %>%
  ggroc() +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") +
  theme_bw() +
  coord_fixed() +
  ggtitle("Union")

roc_mtl <- report_mtl %>% 
  roc(response = "truth", predictor = "GOF") %>%
  ggroc() +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") +
  theme_bw() +
  coord_fixed() +
  ggtitle("MTL")

ggarrange(roc_dirac, roc_union, roc_mtl, roc_mkl)

# PRC
prc_dirac <- report_dirac %>% 
  pr_curve(truth = as.factor(truth), GOF) %>%
  ggplot(aes(x = recall, y = precision)) +
  geom_path() +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") +
  theme_bw() +
  coord_fixed() +
  ggtitle("Dirac")

prc_union <- report_union %>% 
  pr_curve(truth = as.factor(truth), GOF) %>%
  ggplot(aes(x = recall, y = precision)) +
  geom_path() +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") +
  theme_bw() +
  coord_fixed() +
  ggtitle("Union")

prc_mtl <- report_mtl %>% 
  pr_curve(truth = as.factor(truth), GOF) %>%
  ggplot(aes(x = recall, y = precision)) +
  geom_path() +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") +
  theme_bw() +
  coord_fixed() +
  ggtitle("MTL")

prc_mkl <- report_mkl %>% 
  pr_curve(truth = as.factor(truth), GOF) %>%
  ggplot(aes(x = recall, y = precision)) +
  geom_path() +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") +
  theme_bw() +
  coord_fixed() +
  ggtitle("MKL")

ggarrange(prc_dirac, prc_union, prc_mtl, prc_mkl)


