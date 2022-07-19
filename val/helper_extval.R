# helper ext val
# takes a curated list of variants and tool performance for external validation

# pkg
library(librarian)
librarian::shelf(tidyverse,
                 readxl,
                 stringr,
                 hablar)

# data
df <- read_excel("val/ext val 13072022 copy.xlsx")
y <- df$effect

# clean
funcion <- str_split(df$pred_funcion, " ", 2, simplify = T) %>%
  as_tibble(.name_repair = ~c("yhat_funcion", "prob_funcion")) %>%
  retype() %>%
  mutate(y = y)

mtl <- str_split(df$pred_mtl, " ", 2, simplify = T) %>%
  as_tibble(.name_repair = ~c("yhat_mtl", "prob_mtl")) %>%
  retype() %>%
  mutate(y = y)

mtmkl <- str_split(df$pred_mtmkl, " ", 2, simplify = T) %>%
  as_tibble(.name_repair = ~c("yhat_mtmkl", "prob_mtmkl")) %>%
  retype() %>%
  mutate(y = y)

# fix prob (must be prob GOF)
funcion <- funcion %>% 
  mutate(prob_funcion = case_when(yhat_funcion == "LOF" ~ 1-prob_funcion, 
                                  TRUE ~ as.numeric(prob_funcion)))

mtl <- mtl %>% 
  mutate(prob_mtl = case_when(yhat_mtl == "LOF" ~ 1-prob_mtl, 
                                  TRUE ~ as.numeric(prob_mtl)))

mtmkl <- mtmkl %>% 
  mutate(prob_mtmkl = case_when(yhat_mtmkl == "LOF" ~ 1-prob_mtmkl, 
                                  TRUE ~ as.numeric(prob_mtmkl)))

# calc metrics
class_metrics <- metric_set(accuracy, kap, mcc, f_meas, roc_auc, pr_auc)

funcion %>% class_metrics(truth = as.factor(y), prob_funcion, estimate = as.factor(yhat_funcion))
mtl %>% class_metrics(truth = as.factor(y), prob_mtl, estimate = as.factor(yhat_mtl))
mtmkl %>% class_metrics(truth = as.factor(y), prob_mtmkl, estimate = as.factor(yhat_mtmkl))
# WIP probabilistic metrics
