# SCN
# Functional variant prediction for voltage-gated sodium channels
# mtmkl.R

### 05/APR/2022 extremely WIP

# expected
# M_mkl # list of kernel matrices [list(Kt, hpo, M)], each matrix is a view
# y_mkl # vector of labels
# t_vec # vector of task membership

### unregularized group-level weights
constructGroupMKL <- function(x){
  
  # data features
  names_tasks <- unique(t_vec)
  n_tasks <- length(t_vec)
  
  # generate a named list of task indices
  indices_tasks <- vector(mode = "list", length = length(names_tasks))
  for (i in 1:length(names_tasks)) {
    indices_tasks[[i]] <- which(t_vec == names_tasks[[i]])
  }
  names(indices_tasks) <- names_tasks
  
  # generate a list of list with the views for each task
  m_tasks <- list()
  for (i in 1:length(names_tasks)) {
    m_tasks[[i]] <- list(hpo[indices_tasks[[i]], indices_tasks[[i]]],
                         M[indices_tasks[[i]], indices_tasks[[i]]])
  }
  
  # run wrapper MKL for each task
  w_tasks <- list()
  y_d <- vector()
  for (i in 1:length(names_tasks)) {
    y_d <- y_mkl[indices_tasks[[i]]]
    
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
  M <- as.matrix(M)
  
  return(M)
}

### Kandemir method for MTMKL
# cf https://github.com/mitmedialab/PersonalizedMultitaskLearning
constructMTMKL <- function(x, V, C){
  
  # data features
  n_tasks <- length(t_vec)
  n_views <- length(M_mkl)
  
  # parameters
  V <- 0.1 # weight on regularization, Kandemir et al. recommends testing a range from 10^-4 to 10^4
  C <- 100 # cost parameter for SVM
  
  # TODO
  
}
