# SCN
# Functional variant prediction for voltage-gated sodium channels
# mtmkl.R

#' @params M_mkl list of kernel matrices [list(Kt, hpo, M)], where each matrix is a view
#' @params y_mkl vector of labels
#' @params t_vec vector of task membership

x <- M_mkl[2:3]
V <- 0.1 
C <- 100

# reticulate approach?
library(reticulate)
# source_python("LSSVM.py")
# source_python("helperFuncs.py")
# source_python("generic_wrapper.py")
# source_python("wip/MTMKL.py")

### Kandemir method for MTMKL
# cf https://github.com/mitmedialab/PersonalizedMultitaskLearning
constructMTMKL <- function(x, 
                           V, # weight on regularization, Kandemir et al. recommends testing a range from 10^-4 to 10^4
                           C, # cost parameter for SVM
                           debug = TRUE, 
                           regularizer = "none", 
                           max_iter = 100, tol = 0.001){

  # parameters
  self.V <- V
  self.C <- C
  
  self.eta <- array(data = (1/n_views), dim = c(n_tasks, n_views))
  self.last_eta <- self.eta
  
  self.n_instances <- dim(x[[1]])[1]
  self.n_tasks <- length(unique(t_vec))
  self.n_views <- length(x)
  
  if(debug == TRUE) {print(paste0("Initialized with ", self.n_tasks, " tasks and ", self.n_views, " views."))}

  # design choice: L2 regularizer function
  self.regularizer_func <- function(self){
    self.V*sum(V*dist(self.eta, method = "euclidean", diag = TRUE, upper = TRUE))
  }
  
  self.regularizing_grad <- function(self, eta_mat, v, task_index){
    2*v*sum(eta_mat[task_index,]-eta_mat) 
  }

  # TODO
  
}
