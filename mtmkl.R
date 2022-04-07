# SCN
# Functional variant prediction for voltage-gated sodium channels
# mtmkl.R

### 05/APR/2022 extremely WIP

# expected
# M_mkl # list of kernel matrices [list(Kt, hpo, M)], each matrix is a view
# y_mkl # vector of labels
# t_vec # vector of task membership

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
