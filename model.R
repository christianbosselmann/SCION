# SCN
# Functional variant prediction for voltage-gated sodium channels
# model.R

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
                 klic,
                 progress,
                 igraph,
                 dendextend,
                 quiet = TRUE)

#' @params seed random seed
#' @params k inner and outer folds for nested cross validation
#' @params cost_vec vector of cost values to tune over
#' @params class_weight SVM weights for imbalanced classes: "uniform" for equal weight distribution, "inverse" for weights inversely proportional to class frequency
#' @params kernel choice of kernel method: "mkl" for MTMKL-SVM, "mtl" for MTL-SVM, "dirac" for Dirac kernel SVM, "union" for union SVM
#' @params mkl_method choice of MKL method for kernel weights, "semkl" for SEMKL, "simple" for simpleMKL, "uniform" for no kernel weights
#' @params mkl_cost penalty for MKL kernel prioritization (only applies to SEMKL, SimpleMKL)
#' @return saved timestamped objects of parameters, metrics and raw predictions
#' example
#' seed <- 42
#' k <- 10
#' cost_vec <- 2 ^ seq(-5, 5, by = 1)
#' class_weight <- "uniform"
#' kernel <- "mkl"
#' mkl_method <- "uniform"
#' mkl_cost <- 1

# get helper functions
source("func.R")
if (kernel == "rmtl") {source("func_RMTL.R")}
if (kernel == "mkl") {source("func_MKL.R")}

# set random seed
addTaskCallback(function(...) {set.seed(seed);TRUE})

# data features
data <- read_csv("data/dat_prep.csv", col_types = cols())
t_vec <- data$gene 
data$y <- as.factor(data$y)
y <- data$y 

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
  hpo <- kernelPreparation(hpo)
  
  Kt <- readRDS("mat/taskmatrix.rds")
  Kt <- kernelPreparation(Kt)
  
  Km <- readRDS("mat/kernelmatrices_instance.rds")
  Km <- lapply(Km, kernelPreparation)
  
  mkl_weights <- list()
  
  G <- read_csv("mat/distancematrix.csv", col_types = cols()) # graph for RMKL
}
if (kernel == "mtl") {
  Km <- readRDS("mat/kernelmatrices_mtl.rds")
  Km <- unlist(Km, recursive = FALSE)
}
if (kernel == "dirac") {
  Km <- readRDS("mat/kernelmatrices_dirac.rds")
}
if (kernel == "union" || kernel == "rmtl") {
  Km <- readRDS("mat/kernelmatrices_union.rds")
  G <- read_csv("mat/distancematrix.csv", col_types = cols())
}

# set up loop objects
report <- list()
mat_precomp <- list()
mat_splits <- list()
if (mkl_method == "block") {mat_mkl <- list()}

# progress bar
pb <- progress_bar$new(
  format = "(:spin) [:bar] :percent",
  total = length(Km), clear = FALSE, width = 60)

# iterate through the available kernel matrices
for (m in 1:length(Km)) {
  
  # tick progress bar
  pb$tick()
  
  # set up loop objects
  report[[m]] <- data.frame()
  M <- as.matrix(Km[[m]])
  
  if (kernel == "mkl") {M_mkl <- list(Kt, hpo, M)}
  
  # set up nested cv
  cv <- nested_cv(M, outside = vfold_cv(v = k), 
                  inside = vfold_cv(v = k))
  
  # set up cost wrapper functions
  svm_metric <- function(object, cost = 1) {
    i <- object$in_id
    
    if (kernel == "rmtl") {
      d <- constructRMTL(mat = M[i,i], y = y[i], t = t_vec[i], G = G)
      
      M <- applyRMTL(mat = M, t = t_vec, d = d)
      
      assign("M", M, envir = .GlobalEnv)
    }
    if (kernel == "mkl") {
      M_mkl_train <- lapply(M_mkl, function(x) x[i,i])
      if (mkl_method == "uniform") {
        gamma <- rep(1/length(M_mkl), length(M_mkl))
        M <- Reduce(`+`,Map(`*`, gamma, M_mkl))
        assign("M", M, envir = .GlobalEnv)
        mkl_weights[[m]] <<- gamma
      }
      if (mkl_method == "semkl") {
        tryCatch(expr = {
          gamma <- SEMKL.classification(k = M_mkl_train,
                                        outcome = y_mkl[i],
                                        penalty = mkl_cost)$gamma
          M <- Reduce(`+`,Map(`*`, gamma, M_mkl))
          assign("M", M, envir = .GlobalEnv)
        },
        error = function(x) {
          gamma <- rep(1/length(M_mkl_train), length(M_mkl_train))
          M <- Reduce(`+`,Map(`*`, gamma, M_mkl))
          assign("M", M, envir = .GlobalEnv)
        })
      }
      if (mkl_method == "simple") {
        tryCatch(expr = {
          gamma <- SimpleMKL.classification(k = M_mkl_train,
                                            outcome = y_mkl[i],
                                            penalty = mkl_cost)$gamma
          M <- Reduce(`+`,Map(`*`, gamma, M_mkl))
          assign("M", M, envir = .GlobalEnv)
        },
        error = function(x) {
          gamma <- rep(1/length(M_mkl_train), length(M_mkl_train))
          M <- Reduce(`+`,Map(`*`, gamma, M_mkl))
          assign("M", M, envir = .GlobalEnv)
        })
      }
      if (mkl_method == "group") {
        gamma <- constructGroupMKL(matrices = M_mkl_train[2:3], # Task-wise similarity doesn't include useful information for this method
                                   label = y_mkl[i],
                                   tasks = t_vec[i])$gamma
        
        M <- applyGroupMKL(matrices = M_mkl[2:3],
                           tasks = t_vec,
                           gamma = gamma)
        
        assign("M", M, envir = .GlobalEnv)
      }
      if (mkl_method == "block") {
        mod_mkl <- constructBlockMKL(matrices = M_mkl_train, 
                                     label = y_mkl[i],
                                     tasks = t_vec[i],
                                     hierarchical = TRUE,
                                     graph = G)
        
        assign("mkl", mod_mkl, envir = .GlobalEnv)
        
        M <- applyBlockMKL(matrices = M_mkl,
                           tasks = t_vec,
                           gamma = mod_mkl$gamma,
                           hierarchical = TRUE,
                           graph = mod_mkl$graph,
                           delta = mod_mkl$delta)
        
        assign("M", M, envir = .GlobalEnv)
      }
    }
    
    M_train <- M[i,i]
    y_train <- y[i]
    M_test <- as.kernelMatrix(M[-i,i])
    y_test <- y[-i]
    
    mod <- e1071::svm(
      x = M_train,
      y = y_train,
      cost = cost,
      probability = TRUE,
      class.weights = weights
    )
    
    holdout_pred <- kernlab::predict(mod, M_test) %>%
      as_tibble() %>%
      rename(y_hat = value) %>%
      cbind(y_test)
    
    # set metric here
    accuracy(holdout_pred, truth = y_test, estimate = y_hat)$.estimate 
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
  
  # execute inner resampling loops
  plan(multisession)
  tuning_cv <- future_map(cv$inner_resamples, summarize_tune_cv,
                          .options = furrr_options(seed = seed),
                          .progress = FALSE)
  
  # tuning, min/max metric here
  pooled_inner <- tuning_cv %>% bind_rows
  best_cost <- function(dat) dat[which.max(dat$mean_metric),]
  
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
  
  # keep best matrix and resampling info
  mat_precomp[[m]] <- M
  mat_splits[[m]] <- cv
  if (mkl_method == "block") {mat_mkl[[m]] <- mkl}
}

# get combinations of hyperparameters and include in report
report_params <- do.call(rbind, report)
report_params <- report_params %>%
  .[which.max(report_params$mean),] %>% # min/max metric here
  mutate(best_matrix = which.max(report_params$mean)) %>%
  mutate(kernel = kernel) %>%
  mutate(weights = class_weight)

if (kernel == "mkl") {
  report_params <- report_params %>%
    mutate(mkl_cost = mkl_cost) %>%
    mutate(mkl_method = mkl_method)
}

print(report_params)

# assess model with best combination
M <- mat_precomp[[report_params$best_matrix]]
cv <- mat_splits[[report_params$best_matrix]]
C <- as.numeric(levels(report_params$best_cost))[report_params$best_cost] # careful: implicit coercion

# initialize loop 
report_raw <- list()

for (i in 1:length(cv$splits)) {
  train_indices <- cv$splits[[i]]$in_id 
  
  model <- e1071::svm(
    x = M[c(train_indices), c(train_indices)],
    y = y[c(train_indices)],
    cost = C,
    probability = TRUE,
    class.weights = weights
  )
  
  test <- as.kernelMatrix(M[-train_indices, train_indices])
  
  # assess fold
  report_raw[[i]] <- data.table()
  report_raw[[i]]$pred <- kernlab::predict(model, test)
  report_raw[[i]]$truth <- y[-train_indices]
  
  # also get class probabilities
  prob <- attr(predict(model, test, probability = TRUE), "probabilities")
  report_raw[[i]] <- cbind(report_raw[[i]], prob)
  
  # ...and distance from hyperplane for histogram of projections
  dist <- attr(predict(model, test, decision.value = TRUE), "decision.values") 
  colnames(dist) <- "dist"
  report_raw[[i]] <- cbind(report_raw[[i]], dist)
  
  # store test indices to later get task-specific performance metrics
  report_raw[[i]]$ind <- which(1:nrow(M) %nin% train_indices)
}

# generate reports
generateReport(report_raw, print = TRUE, export = TRUE)
