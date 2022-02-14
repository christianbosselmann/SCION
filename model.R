# SCN
# Functional variant prediction for voltage-gated sodium channels
# model.R

# @params
seed <- 42 # random seed
k <- 5 # k-fold nested cross validation
cost_vec <- 2 ^ seq(-5, 5, by = 1) # C-SVM hyperparameter vector for grid search
class_weight <- "uniform" # SVM weights for imbalanced classes: "uniform" for equal weight distribution, "inverse" for weights inversely proportional to class frequency
kernel <- "mkl" # choice of kernel matrices: "mkl" for MTMKL-SVM, "mtl" for MTL-SVM, "dirac" for Dirac kernel SVM, "union" for union SVM
mkl_method <- "semkl" # choice of MKL method for kernel weights, "semkl" for SEMKL, "simple" for simpleMKL, "uniform" for no kernel weights
mkl_cost <- 1 # penalty for MKL kernel prioritization (SEMKL, SimpleMKL)

# load packages
library("librarian")
librarian::shelf(tidyverse, 
                 tidymodels,
                 data.table,
                 kernlab,
                 furrr,
                 caret,
                 cowplot,
                 assertthat,
                 Matrix,
                 matrixcalc,
                 quiet = TRUE)

# get helper functions
source("func.R")

# set random seed
set.seed(seed)

# read data
data <- read_csv("data/dat_prep.csv")

t_vec <- data$gene # get task vector
data$gene <- NULL # drop task from features

data$y <- as.factor(data$y)
y <- data$y # get label vector

# set up class weights
if (class_weight == "uniform") {
  weights <- c("GOF" = 1,
               "LOF" = 1)
}

if (class_weight == "inverse") {
  weights <- c("GOF" = length(data$y) / sum(data$y == "GOF"),
               "LOF" = length(data$y) / sum(data$y == "LOF"))
}

# load kernel matrices, maintain order, and set up loop objects depending on method
if (kernel == "mkl") {
  y_mkl <- as.numeric(y) # another label vector for MKL
  y_mkl[y_mkl == 1] <- 1
  y_mkl[y_mkl == 2] <- -1
  
  load("mat/hpomatrix.RData")
  hpo <- round(hpo, 10)
  hpo <- normalize.kernel(hpo, method = "sqrt")
  
  Kt <- readRDS("mat/taskmatrix.rds")
  Kt <- round(Kt, 10)
  Kt <- normalize.kernel(Kt, method = "sqrt")
  
  Km <- readRDS("mat/kernelmatrices_instance.rds")
  Km <- lapply(Km, round, 10)
  Km <- lapply(Km, normalize.kernel, method = "sqrt")
}

if (kernel == "mtl") {Km <- readRDS("mat/kernelmatrices_mtl.rds")}
if (kernel == "dirac") {Km <- readRDS("mat/kernelmatrices_dirac.rds")}
if (kernel == "union") {Km <- readRDS("mat/kernelmatrices_union.rds")}

if (!(kernel == "mkl")) {Km <- unlist(Km, recursive = FALSE)}
if (kernel == "mkl") {mkl_weights <- list()}

report <- list()
mat_precomp <- list()

# iterate through the available kernel matrices
for (m in 1:length(Km)) {
  report[[m]] <- data.frame()
  
  M <- as.matrix(Km[[m]])

  ### MKL module
  if (kernel == "mkl") {
    M_mkl <- list(Kt, hpo, M)
    
    if (mkl_method == "simple") {
      mod_mkl <- SimpleMKL.classification(k = M_mkl, 
                                          outcome = y_mkl,
                                          penalty = mkl_cost)
    }
    
    if (mkl_method == "semkl") {
      mod_mkl <- SEMKL.classification(k = M_mkl, 
                                      outcome = y_mkl,
                                      penalty = mkl_cost)
    }
    
    if (mkl_method == "uniform") {
      mod_mkl <- NULL
      mod_mkl$gamma <- c(1/3,1/3,1/3)
    }
    
    # finalize MKL matrix by taking the weighted mean
    M <- Reduce(`+`,Map(`*`, mod_mkl$gamma, M_mkl))
    
    # store weights
    mkl_weights[[m]] <- mod_mkl$gamma
  }
  ###
  
  # set up nested cv
  cv <- nested_cv(M, outside = vfold_cv(v = k), 
                  inside = vfold_cv(v = k))
  
  # set up cost wrapper functions
  svm_metric <- function(object, cost = 1) {
    y_test <- y[-object$in_id]
    
    mod <- e1071::svm(
      x = M[c(object$in_id), c(object$in_id)],
      y = y[c(object$in_id)],
      cost = cost,
      probability = TRUE,
      class.weights = weights
    )
    
    test <- as.kernelMatrix(M[-object$in_id, object$in_id])
    
    holdout_pred <- predict(mod, test) %>%
      as_tibble() %>%
      rename(y_hat = value) %>%
      cbind(y_test)
    
    accuracy(holdout_pred, truth = y_test, estimate = y_hat)$.estimate # set metric here
  }
  
  # parameterize
  metric_wrapper <- function(cost, object) svm_metric(object, cost)
  
  # purrr over cost and add new column 'metric' to store results
  tune_over_cost <- function(object) {
    tibble(cost = cost_vec) %>% 
      mutate(metric = map_dbl(cost, metric_wrapper, object = object))
  }
  
  # this will be called over the set of outer cv splits
  summarize_tune_cv <- function(object) {
    map_df(object$splits, tune_over_cost) %>%
      group_by(cost) %>%
      summarize(mean_metric = mean(metric, na.rm = TRUE),
                n = length(metric),
                .groups = "drop")
  }
  
  # execute inner resampling loops in parallel via furrr
  plan(multisession)
  tuning_cv <- future_map(cv$inner_resamples, summarize_tune_cv, 
                          .options = furrr_options(seed = seed)) 
  
  # tuning
  pooled_inner <- tuning_cv %>% bind_rows
  best_cost <- function(dat) dat[which.max(dat$mean_metric),] # min/max metric here
  
  # best parameter estimate from outer resampling
  cost_vals <- tuning_cv %>% 
    map_df(best_cost) %>% 
    select(cost)
  
  cv <- bind_cols(cv, cost_vals) %>% 
    mutate(cost = factor(cost, levels = paste(cost_vec)))
  
  # compute the outer resampling results for each of the splits using the tuning val
  cv <- cv %>% 
    mutate(metric = map2_dbl(splits, cost, svm_metric))
  
  # generate report
  report[[m]] <- tibble(mean = mean(cv$metric),
                        sd = sd(cv$metric),
                        best_cost = cv$cost[which.max(cv$metric)]) # min/max metric here
  
  # keep best matrix
  mat_precomp[[m]] <- M
}

# get best combinations of hyperparameters
report_params <- do.call(rbind, report)
report_params <- report_params %>%
  .[which.max(report_params$mean),] %>% # min/max metric here
  mutate(best_matrix = which.max(report_params$mean)) %>%
  mutate(kernel = kernel) %>%
  mutate(weights = class_weight)

print(report_params)

# assess model with best combination
M2 <- mat_precomp[[report_params$best_matrix]]
C <- as.numeric(levels(report_params$best_cost))[report_params$best_cost] # careful: implicit coercion

# initialize loop object
report_raw <- list()

for (i in 1:length(cv$splits)) {
  # training
  train_indices <- cv$splits[[i]]$in_id 
  
  model <- e1071::svm(
    x = M2[c(train_indices), c(train_indices)],
    y = y[c(train_indices)],
    cost = C,
    probability = TRUE,
    class.weights = weights
  )
  
  test <- as.kernelMatrix(M2[-train_indices, train_indices])
  
  # assess fold
  report_raw[[i]] <- data.table()
  report_raw[[i]]$pred <- kernlab::predict(model, test)
  report_raw[[i]]$truth <- y[-train_indices]
  
  # also get class probabilities
  prob <- attr(predict(model, test, probability = TRUE), "probabilities")
  report_raw[[i]] <- cbind(report_raw[[i]], prob)
  
  # store test indices to later get task-specific performance metrics
  report_raw[[i]]$ind <- which(1:nrow(M2) %nin% train_indices) # very ugly hack to get test indices
}

# keep raw predictions and get metrics
class_metrics <- metric_set(accuracy, kap, mcc, f_meas, roc_auc, pr_auc)

report_predictions <- report_raw %>%
  rbindlist(., idcol = "fold", use.names = TRUE) %>%
  group_by(fold) %>%
  class_metrics(truth = truth, GOF, estimate = pred) %>%
  group_by(.metric) %>%
  summarise(mean = mean(.estimate, na.rm = TRUE), sd = sd(.estimate, na.rm = TRUE))

ind_list <- t_vec %>% # list of task indices for report
  as_tibble() %>%
  tibble::rowid_to_column("ind") %>%
  rename(gene = value)

report_final <- report_raw %>% 
  rbindlist(., idcol = "fold", use.names = TRUE) %>%
  inner_join(ind_list, by = "ind")

print(report_predictions)

# export model hyperparameters, raw predictions and model performance
write_csv(report_params, paste0('out/report_params_', Sys.time(), '.csv'))
write_csv(report_final, paste0('out/report_preds_', Sys.time(), '.csv'))
write_csv(report_predictions, paste0('out/report_metrics_', Sys.time(), '.csv'))
