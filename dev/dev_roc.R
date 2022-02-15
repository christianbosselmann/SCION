# SCN dev
# this script takes "report_final.csv" files, i.e. the raw prediction from nested k-fold cv, and prepares a grid of ROC plots
# just a code snippet, remove before publication

library(librarian)
librarian::shelf(ggplot2,
                 ggpubr,
                 grid)

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


