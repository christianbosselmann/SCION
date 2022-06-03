# SCN
# Functional variant prediction for voltage-gated sodium channels
# some helper functions

### negate in
`%nin%` = Negate(`%in%`)

### label conversion for MKL functions
makeLabelNumeric <- function(y){
  if(!is.numeric(y)) {
    y <- as.numeric(y) # implicit boolean conversion
    y[y == 1] <- 1 
    y[y == 2] <- -1
    return(y)
  }
  if(is.numeric(y)) {
    return(y)
  }
}

### generating reports from predictions
generateReport <- function(x, print = FALSE, export = FALSE){
  class_metrics <- metric_set(accuracy, kap, mcc, f_meas, roc_auc, pr_auc)
  
  report_predictions <- report_raw %>%
    rbindlist(., idcol = "fold", use.names = TRUE) %>%
    group_by(fold) %>%
    class_metrics(truth = truth, GOF, estimate = pred) %>%
    group_by(.metric) %>%
    summarise(mean = mean(.estimate, na.rm = TRUE), sd = sd(.estimate, na.rm = TRUE))
  
  ind_list <- t_vec %>% 
    as_tibble() %>%
    tibble::rowid_to_column("ind") %>%
    rename(gene = value)
  
  report_final <- report_raw %>% 
    rbindlist(., idcol = "fold", use.names = TRUE) %>%
    inner_join(ind_list, by = "ind")
  
  if(print == TRUE){
    print(report_predictions)
  }
  
  if(export == TRUE){
    write_csv(report_params, paste0('out/report_params_', Sys.time(), '.csv'))
    write_csv(report_final, paste0('out/report_preds_', Sys.time(), '.csv'))
    write_csv(report_predictions, paste0('out/report_metrics_', Sys.time(), '.csv'))
  }
  
}

### normalize kernel matrix
# cf. Kernel Methods for Pattern Analysis, Algorithm 5.1
#' @param K kernel matrix to be normalized
#' @return normalized kernel matrix
kernelNormalisation <- function(K){
  # min_eigenvalue <- min(eigen(K)$values)
  # K <- sqrt(diag(K) + min_eigenvalue)
  D <- diag(1/sqrt(diag(K)))
  K <- D %*% K %*% D
  return(K)
}

### center in feature space
# cf. Kernel Methods for Pattern Analysis, Algorithm 5.3
#' @param K kernel matrix to be centered
#' @return kernel matrix centered in feature space
kernelCentering <- function(K){
  ell <- dim(K)[1]
  D <- colSums(K)/ell # row vector storing the column averages of K
  E <- sum(D)/ell # average of all the entries of K
  J <- matrix(1, ell, 1) %*% D
  Jt <- Conj(t(J)) # complex conjugate transpose of J
  K <- K - J - Jt + E * matrix(1, ell, ell)
}

### check if matrix is symmetric and psd
#' @param K kernel matrix to be checked
#' @return prints matrix properties to console
#' adapted from base and matrixcalc
kernelCheck <- function(K, tol = 1e-08){
  if(!isSymmetric(K)){stop("Argument is not a symmetric matrix.")}
  if(isSymmetric(K)){print("Argument is a symmetric matrix.")}
  
  ev <- eigen(K)[[1]]
  n <- nrow(K)
  
  for (i in 1:n) {
    if (abs(ev[i]) < tol) {
      ev[i] <- 0
    }
  }
  if (any(ev < 0)) {
    stop("Argument is not a psd matrix.")
  }
  print("Argument is a psd matrix.")
}

### standard kernel preprocessing
#' @param K kernel matrix to be normalized and centered
#' @return kernel matrix 
kernelPreparation <- function(K){
  K <- kernelNormalisation(K)
  K <- kernelCentering(K)
  K <- round(K, 10)
  return(K)
}

#' Extract step item
#' Returns extracted step item from prepped recipe.
#' attr. Giovanni Colitti
#' @param recipe Prepped recipe object.
#' @param step Step from prepped recipe.
#' @param item Item from prepped recipe.
#' @param enframe Should the step item be enframed?
#' @export
extract_step_item <- function(recipe, step, item, enframe = TRUE) {
  d <- recipe$steps[[which(purrr::map_chr(recipe$steps, ~ class(.)[1]) == step)]][[item]]
  if (enframe) {
    tibble::enframe(d) %>% tidyr::spread(key = 1, value = 2)
  } else {
    d
  }
}

#' Unnormalize variable
#' Unnormalizes variable using standard deviation and mean from a recipe object. See \code{?recipes}.
#' attr. Giovanni Colitti
#' @param x Numeric vector to normalize.
#' @param rec Recipe object.
#' @param var Variable name in the recipe object.
#' @export
unNormalize <- function(x, rec, var) {
  var_sd <- extract_step_item(rec, "step_normalize", "sds") %>% dplyr::pull(var)
  var_mean <- extract_step_item(rec, "step_normalize", "means") %>% dplyr::pull(var)
  
  (x * var_sd) + var_mean
}
