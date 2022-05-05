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
generateReport <- function(x, print, export){
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

### unregularized group-level weights (group MKL, or Dirac MKL)
#' @param matrices list of Kernel matrices
#' @param tasks vector of task membership
#' @param label MKL label vector
#' @return mod_mkl containg M (training kernel matrix) and gamma (weight vector)
constructGroupMKL <- function(matrices, label, tasks){
  # tasks <- t_vec
  # matrices <- M_mkl
  # gamma <- c(0.3, 0.3, 0.3)
  # label <- y_mkl
  
  mod_mkl <- list()
  
  # data features
  names_tasks <- unique(tasks)
  
  # generate a named list of task indices
  indices_tasks <- vector(mode = "list", length = length(names_tasks))
  for (i in 1:length(names_tasks)) {
    indices_tasks[[i]] <- which(tasks == names_tasks[[i]])
  }
  names(indices_tasks) <- names_tasks
  
  # generate a list of list with the views for each task
  m_tasks <- list()
  for (i in 1:length(names_tasks)) {
    m_tasks[[i]] <- lapply(matrices, function(x) x[indices_tasks[[i]], indices_tasks[[i]]])
  }
  
  # run wrapper MKL for each task
  w_tasks <- list()
  y_d <- vector()
  for (i in 1:length(names_tasks)) {
    y_d <- label[indices_tasks[[i]]] 
    
    # small tasks with similar observations may lead to computationally singular systems
    # if this occurs, return the uniformly weighted kernel matrix
    tryCatch(
      expr = {
        w_tasks[[i]] <- SEMKL.classification(k = m_tasks[[i]],
                                             outcome = y_d,
                                             penalty = mkl_cost)$gamma
      },
      error = function(x) {
        w_tasks[[i]] <<- rep(1/length(m_tasks[[i]]), length(m_tasks[[i]]))
      })
  }
  
  # apply weights and store in the original list of matrices
  for (i in 1:length(names_tasks)) {
    m_tasks[[i]] <- Reduce(`+`,Map(`*`, w_tasks[[i]], m_tasks[[i]]))
  }
  
  # merge into a diagonal block matrix
  M <- bdiag(m_tasks)
  M <- as.data.frame(as.matrix(M))
  
  # reorder matrix to match original indices
  indices_tasks <- as.numeric(unlist(indices_tasks))
  M <- cbind(M, indices_tasks)
  M <- arrange(M, indices_tasks) %>%
    select(-indices_tasks)
  M <- as.matrix(M)
  
  # store objects for output
  mod_mkl$M <- M
  mod_mkl$gamma <- w_tasks
  
  return(mod_mkl)
}

#' @param matrices list of list of Kernel matrices, where each sublist is a view
#' @param tasks vector of task membership
#' @param gamma list of view weight vectors from constructGroupMKL function
#' @return mod_mkl containg M (training kernel matrix) and gamma (weight vector)
applyGroupMKL <- function(matrices, tasks, gamma) {
  
  # data features
  names_tasks <- unique(tasks)
  
  # generate a named list of task indices
  indices_tasks <- vector(mode = "list", length = length(names_tasks))
  for (i in 1:length(names_tasks)) {
    indices_tasks[[i]] <- which(tasks == names_tasks[[i]])
  }
  names(indices_tasks) <- names_tasks
  
  # generate a list of list with the views for each task
  m_tasks <- list()
  for (i in 1:length(names_tasks)) {
    m_tasks[[i]] <- lapply(matrices, function(x) x[indices_tasks[[i]], indices_tasks[[i]]])
  }
  
  # apply weights
  m_tasks <- lapply(seq_along(m_tasks), function (x) Reduce(`+`,Map(`*`, gamma[[x]], m_tasks[[x]])))
  
  # merge into a diagonal block matrix
  M <- bdiag(m_tasks)
  M <- as.data.frame(as.matrix(M))
  
  # reorder matrix to match original indices and return
  indices_tasks <- as.numeric(unlist(indices_tasks))
  M <- cbind(M, indices_tasks)
  M <- arrange(M, indices_tasks) %>%
    select(-indices_tasks)
  M <- as.matrix(M)
  
  return(M)
}

### blockwise MTMKL
#' @param matrices list of Kernel matrices, where the first element is Kt (task-wise similarity)
#' @param label label vector of training observations
#' @param tasks task membership vector t_vec of training observations
#' @param graph network graph or similarity/distance matrix of dim TxT
#' @return M training MKL matrix
#' cf. https://www.cmpe.boun.edu.tr/~ethem/files/papers/gonen_icml08.pdf
#' cf. DOI 10.1186/1471-2105-11-S8-S5
#' this method decomposes the problem into t*t MKL learning problems, where weights are learned for each combination of tasks (leaves of the hierarchy) which reduces to a quasi-conformal transformation (Amari & Wu, 1998)
constructBlockMKL <- function(matrices, label, tasks, graph = NULL){
  
  # dependencies
  requireNamespace("Matrix", quietly = TRUE)
  requireNamespace("tidyverse", quietly = TRUE)
  # requireNamespace("parallel", quietly = TRUE)
  
  # split into task-wise similarity matrix and instance-level matrices to weigh
  Kt <- matrices[[1]]
  matrices <- matrices[2:3]
  
  # data features
  dim <- dim(matrices[[1]])[1]
  names_tasks <- unique(tasks)
  y <- makeLabelNumeric(label)
  
  # generate a named list of task indices
  indices_tasks <- vector(mode = "list", length = length(names_tasks))
  for (i in 1:length(names_tasks)) {indices_tasks[[i]] <- which(tasks == names_tasks[[i]])}
  names(indices_tasks) <- names_tasks
  
  # get new task indices
  indices_tasks_new <- list()
  for (i in 1:length(names_tasks)) {
    l <- length(indices_tasks[[i]]) # length of index vector
    k <- length(unlist(indices_tasks[(1:i)-1])) # length of all previous tasks index vectors
    t <- c(k+(1:l)) # task block coordinates
    indices_tasks_new[[i]] <- t
  }
  names(indices_tasks_new) <- names_tasks
  
  list_2tasks <- expand.grid(names_tasks, names_tasks, stringsAsFactors = FALSE) %>%
    as.data.frame() %>%
    split(., seq(nrow(.)), drop = TRUE)
  
  m_2tasks <- list() # list of combined task matrices
  wt <- list() # list of weight vectors

  for (i in 1:length(list_2tasks)) {
    i1 <- indices_tasks[[list_2tasks[[i]][[1]]]] # indices of first task to grab from
    i2 <- indices_tasks[[list_2tasks[[i]][[2]]]] # indices of second task to grab from

    m <- lapply(matrices, function(x) x[c(i1,i2), c(i1,i2)])
    y_d <- y[c(i1,i2)]

    # small tasks with similar observations may lead to computationally singular systems
    # if this occurs, return the uniformly weighted kernel matrix
    tryCatch(expr = {
      w <- SEMKL.classification(k = m,
                                outcome = y_d,
                                penalty = mkl_cost)$gamma
    },
    error = function(x) {
      w <<- rep(1/length(m), length(m))
    })

    # store weights of task combinations
    wt[[i]] <- w

    # merge view matrices by weights
    m_2tasks[[i]] <- Reduce(`+`,Map(`*`, w, m))
  }

  # name stored weight vector
  wt <- do.call(rbind.data.frame, wt)
  colnames(wt) <- c("hpo", "M")
  wt$t1 <- unlist(lapply(lapply(list_2tasks, function (x) x[1,]), function (x) paste(x[,1])))
  wt$t2 <- unlist(lapply(lapply(list_2tasks, function (x) x[1,]), function (x) paste(x[,2])))
  # assign("wt", wt, envir = .GlobalEnv)
  
  # TODO fix: "x replacement has 81 rows, data has 71"
  # # parallel version of the above for loop, speedup 22s to 6.8s
  # p <- list()
  # p <- mclapply(list_2tasks, parallelBlockMKL, 
  #               ind = indices_tasks, mat = matrices, y = y, mkl_cost = 1,
  #               mc.cores = detectCores())
  # 
  # # extract and name stored weight vector
  # wt <- lapply(p, function(x) rbind(x[1]$w)) %>%
  #   do.call(rbind,.) %>%
  #   as.data.frame()
  # colnames(wt) <- c("hpo", "M")
  # wt$t1 <- NA
  # wt$t1 <- unlist(lapply(lapply(list_2tasks, function (x) x[1,]), function (x) paste(x[,1])))
  # wt$t2 <- NA
  # wt$t2 <- unlist(lapply(lapply(list_2tasks, function (x) x[1,]), function (x) paste(x[,2])))
  # assign("wt", wt, envir = .GlobalEnv)
  # 
  # # extract and store list of task combinations per matrix
  # m_2tasks <- lapply(p, function(x) x[2]$m2)
  
  # preallocate matrix
  M <- matrix(data = 0, nrow = dim, ncol = dim)
  
  # populate matrix with taskwise weighted MKL matrices
  # [row,col]
  for (i in 1:length(list_2tasks)){
    j1 <- indices_tasks_new[[list_2tasks[[i]][[1]]]] # indices of first task to store later
    j2 <- indices_tasks_new[[list_2tasks[[i]][[2]]]] # indices of second task to store later
    
    m <- m_2tasks[[i]]
    
    k1 <- 1:length(j1)
    k2 <- 1:length(j2)
    mul <- m[k2,k2]
    mur <- m[k2,k1]
    mll <- m[k1,k2]
    mlr <- m[k1,k1]
    
    M[j2,j2] <- M[j2,j2] + mul
    M[j2,j1] <- M[j2,j1] + mur
    M[j1,j2] <- M[j1,j2] + mll
    M[j1,j1] <- M[j1,j1] + mlr
  }
  
  # symmetric matrix reorder 
  M <- M[order(unlist(indices_tasks, use.names = FALSE)),]
  M <- M[,order(unlist(indices_tasks, use.names = FALSE))]
  
  # get product kernel matrix with task similarity matrix for MTMKL
  M <- M * Kt
  
  if(!is.null(graph)){
    ### TODO hierarchical decomposition
    # # set up hierarchy
    # graph <- read_csv("mat/distancematrix.csv", col_types = cols())
    # tmp <- as.dist(graph)
    # tmp <- hclust(tmp, method = "complete")
    # plot(tmp)
  }
  
  # store output
  d <- list()
  d$M <- M # training kernel matrix
  d$gamma <- wt # task-wise weights
  
  return(d)
} 

#' @param matrices list of Kernel matrices, where the first element is Kt (task-wise similarity)
#' @param label label vector 
#' @param tasks vector t_vec of task membership 
#' @param gamma task-wise weights learned in constructBlockMKL fn
#' @return M MKL for training/testing
#' this function applies weights learned from the training indices and yields a complete MKL matrix for train/test splits
applyBlockMKL <- function(matrices, tasks, gamma) {
  # debug
  # matrices <- M_mkl
  # tasks <- t_vec
  # gamma <- gamma
  
  # split into task-wise similarity matrix and instance-level matrices to weigh
  Kt <- matrices[[1]]
  matrices <- matrices[2:3]
  
  # data features
  dim <- dim(matrices[[1]])[1]
  names_tasks <- unique(tasks)
  
  # generate a named list of task indices
  indices_tasks <- vector(mode = "list", length = length(names_tasks))
  for (i in 1:length(names_tasks)) {indices_tasks[[i]] <- which(tasks == names_tasks[[i]])}
  names(indices_tasks) <- names_tasks
  
  # get new task indices
  indices_tasks_new <- list()
  for (i in 1:length(names_tasks)) {
    l <- length(indices_tasks[[i]]) # length of index vector
    k <- length(unlist(indices_tasks[(1:i)-1])) # length of all previous tasks index vectors
    t <- c(k+(1:l)) # task block coordinates
    indices_tasks_new[[i]] <- t
  }
  names(indices_tasks_new) <- names_tasks
  
  list_2tasks <- expand.grid(names_tasks, names_tasks, stringsAsFactors = FALSE) %>%
    as.data.frame() %>%
    split(., seq(nrow(.)), drop = TRUE)
  
  m_2tasks <- list() # list of combined task matrices
  for (i in 1:length(list_2tasks)) {
    i1 <- indices_tasks[[list_2tasks[[i]][[1]]]] # indices of first task to grab from
    i2 <- indices_tasks[[list_2tasks[[i]][[2]]]] # indices of second task to grab from
    m <- lapply(matrices, function(x) x[c(i1,i2), c(i1,i2)])
    
    # select appropriate row from precomputed weight vector of constructBlockMKL fn
    w <- as.numeric(gamma[i,][1:2]) # TODO is this correct?
    
    # merge view matrices by weights
    m_2tasks[[i]] <- Reduce(`+`,Map(`*`, w, m))
  }
  
  # preallocate matrix
  M <- matrix(data = 0, nrow = dim, ncol = dim)
  
  # populate matrix with taskwise weighted MKL matrices
  # [row,col]
  for (i in 1:length(list_2tasks)){
    j1 <- indices_tasks_new[[list_2tasks[[i]][[1]]]] # indices of first task to store later
    j2 <- indices_tasks_new[[list_2tasks[[i]][[2]]]] # indices of second task to store later
    
    m <- m_2tasks[[i]]
    
    k1 <- 1:length(j1)
    k2 <- 1:length(j2)
    mul <- m[k2,k2]
    mur <- m[k2,k1]
    mll <- m[k1,k2]
    mlr <- m[k1,k1]
    
    M[j2,j2] <- M[j2,j2] + mul
    M[j2,j1] <- M[j2,j1] + mur
    M[j1,j2] <- M[j1,j2] + mll
    M[j1,j1] <- M[j1,j1] + mlr
  }
  
  # symmetric matrix reorder 
  M <- M[order(unlist(indices_tasks, use.names = FALSE)),]
  M <- M[,order(unlist(indices_tasks, use.names = FALSE))]
  
  # get product kernel matrix with task similarity matrix for MTMKL
  M <- M * Kt
  
  return(M)
}

# helper function to parallelize block-wise weight learning for constructBlockMKL
parallelBlockMKL <- function(x = list_2tasks, ind = indices_tasks, mat = matrices, y = y, mkl_cost = 1){
  p <- list()
  
  i1 <- ind[[x[[1]]]]
  i2 <- ind[[x[[2]]]]
  
  m <- lapply(mat, function(i) i[c(i1,i2), c(i1,i2)])
  y <- y[c(i1,i2)]
  
  # small tasks with similar observations may lead to computationally singular systems
  # if this occurs, return the uniformly weighted kernel matrix
  tryCatch(expr = {
    p$w <- SEMKL.classification(k = m,
                                outcome = y,
                                penalty = mkl_cost)$gamma
  },
  error = function(x) {
    p$w <<- rep(1/length(m), length(m))
  })
  
  # merge view matrices by weights
  p$m2 <- Reduce(`+`,Map(`*`, p$w, m))
  
  return(p)
}

