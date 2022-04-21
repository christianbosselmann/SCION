# SCN
# Functional variant prediction for voltage-gated sodium channels
# some helper functions

### negate in
`%nin%` = Negate(`%in%`)

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
#' adapted from base and matrixcalc, just for debugging
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

### SimpleMKL, from RMKL pkg
#' Simple MKL 
#'
#' This function conducts Simple MKL for precomputed gramm matrices
#' @param k list of Gramm matrices
#' @param outcome vector of binary outcome -1 and 1
#' @param penalty ppenalty of the smoothness of the resulting desicion rules
#' @param tol change between to iterations is smaller than this, algorithms is considered to have converged
#' @param max.iters maximum number of allowed iteratons
#' @return gamma weight vector for the importnace of each kernel 
#' @return alpha coeffiencents of the dual of MKL
#' @return time total amount of time to train model
#' @return max.iters Numvber of iterations to reach convergence criteria
#' @export
#' @examples
#' library(kernlab)
#' library(caret)
#' library(RMKL)
#' #Load data
#' data(benchmark.data)
#' example.data=benchmark.data[[1]]
#' # Split samples into training and test sets 
#' training.samples=sample(1:dim(example.data)[1],floor(0.7*dim(example.data)[1]),replace=FALSE)
#' # Set up cost parameters and kernels 
#' C=100
#' kernels=rep('radial',3)
#' degree=rep(0,3)
#' scale=rep(0,3)
#' sigma=c(0,2^seq(-3:0))
#' K=kernels.gen(example.data[,1:2], training.samples, kernels, degree, scale, sigma)
#' K.train=K$K.train
#' SimpleMKL.classification(K.train,example.data[training.samples,3], C)

SimpleMKL.classification=function(k,outcome,penalty,tol=10^(-4),max.iters=1000){
  inner=function(Ddag,gammadag,k,J,penalty){
    Jdag=0
    while(Jdag<J&sum(Ddag<0)>0){
      mu=min(which(gammadag==max(gammadag)))
      D=Ddag/sum(abs(Ddag))
      gamma=gammadag
      neg=which(D<0&gamma>0)
      v=neg[which(-gamma[neg]/D[neg]==min(-gamma[neg]/D[neg],na.rm=TRUE))]
      stepsize=-gamma[v]/D[v]
      gammadag=gamma+stepsize*D
      Ddag[mu]=D[mu]-D[v]
      Ddag[v]=0
      kk=Reduce('+',mapply("*", k, gammadag,SIMPLIFY = FALSE))
      h=kk*(outcome%*%t(outcome))
      model=kernlab::ipop(rep(-1,length(outcome)),h,t(outcome),0,rep(0,length(outcome)),rep(penalty,length(outcome)),0)
      Jdag=-1/2*sum(h*kernlab::primal(model)%*%t(kernlab::primal(model)))+sum(kernlab::primal(model))
    }
    return(list('gamma'=gamma,'objval'=Jdag,'direction'=D,'stepsize'=stepsize))
  }
  
  gamma_all=list()
  gamma=rep(1/length(k),length(k))
  epsilon=1
  iters=0
  while(epsilon>tol&iters<max.iters&sum(gamma==1)==0){
    # tic()
    iters=iters+1
    gamma_all[[iters]]=gamma
    kk=Reduce('+',mapply("*", k, gamma,SIMPLIFY = FALSE))
    h=kk*(outcome%*%t(outcome))
    model=kernlab::ipop(rep(-1,length(outcome)),h,t(outcome),0,rep(0,length(outcome)),rep(penalty,length(outcome)),0)
    J=-1/2*sum(h*kernlab::primal(model)%*%t(kernlab::primal(model)))+sum(kernlab::primal(model))
    dJ=sapply(1:length(k), function(a) -1/2*sum(k[[a]]*outcome%*%t(outcome)*kernlab::primal(model)%*%t(kernlab::primal(model))))
    mu=min(which(gamma==max(gamma)))
    gradJ=rep(0,length(k))
    gradJ[-mu]=dJ[-mu]-dJ[mu]
    gradJ[mu]=-sum(gradJ)
    cond1=which(gamma==0&gradJ>0)
    cond2=setdiff(which(gamma>0),mu)
    cond3=mu
    D=rep(0,length(k))
    D[cond1]=0
    D[cond2]=-gradJ[cond2]
    D[cond3]=sum(gradJ[cond2])
    Ddag=D
    innersol=inner(Ddag=D,gammadag = gamma,k=k,J,penalty)
    obj=function(gamma,outcome,penalty,direction,stepsize){
      gammanew=gamma+direction*stepsize
      kk=Reduce('+',mapply("*", k, gammanew,SIMPLIFY = FALSE))
      h=kk*(outcome%*%t(outcome))
      model=kernlab::ipop(rep(-1,length(outcome)),h,t(outcome),0,rep(0,length(outcome)),rep(penalty,length(outcome)),0)
      J=-1/2*sum(h*kernlab::primal(model)%*%t(kernlab::primal(model)))+sum(kernlab::primal(model))
      return(J)
    }
    step_opt=stats::optim(innersol$stepsize/2,obj,outcome=outcome,penalty=penalty,direction=innersol$direction,gamma=innersol$gamma,lower=0,upper=innersol$stepsize,method='L-BFGS-B')
    epsilon=max(gamma-innersol$gamma+step_opt$par*innersol$direction)
    gamma=innersol$gamma+step_opt$par*innersol$direction
    gamma=gamma/sum(gamma)
    #toc(log = TRUE, quiet = TRUE)
  }
  gamma_all[[iters+1]]=gamma
  # log.lst <- tic.log(format = FALSE)
  j=match(kernlab::primal(model)[(kernlab::primal(model)>0)&(kernlab::primal(model)<penalty)][1],kernlab::primal(model))
  b=outcome[j]-sum(kernlab::primal(model)*outcome*kk[,j])
  #tic.clearlog()
  return(list('gamma'=gamma,'iters'=iters,'alpha'=kernlab::primal(model),'b'=b,'gamma_all'=gamma_all))
}

### SEMKL, from RMKL pkg
SEMKL.classification=function(k,outcome,penalty,tol=0.0001,max.iters=1000){
  delta=rep(1,length(k))
  iters=0
  n=length(outcome)
  m=length(k)
  gamma=rep(1/length(k),length(k))
  gamma_all=list()
  #tic()
  while (max(delta)>tol && iters<max.iters){
    iters=iters+1
    gamma_all[[iters]]=gamma
    kg=lapply(1:length(k), function(a) k[[a]]*gamma[a])
    kk=Reduce('+',kg)
    h=kk*(outcome%*%t(outcome))
    model=kernlab::ipop(rep(-1,n),h,t(outcome),0,rep(0,n),rep(penalty,n),0)
    alpha=kernlab::primal(model)
    fnorm=sapply(1:length(k), function(a) sqrt(gamma[a]^2*(alpha*outcome)%*%k[[a]]%*%(alpha*outcome)))
    temp=gamma
    gamma=fnorm/sum(fnorm)
    delta=abs(temp-gamma)
  }
  #toc(log=TRUE,quiet=TRUE)
  #time=tic.log(format=FALSE)[[1]]$toc-tic.log(format=FALSE)[[1]]$tic
  j=match(alpha[(alpha>0)&(alpha<penalty*0.9999)][1],alpha)
  b=outcome[j]-sum(alpha*outcome*kk[,j])
  gamma_all[[iters+1]]=gamma
  results=list("alpha"=alpha,"b"=b,"gamma"=temp,"iters"=iters,'gamma_all'=gamma_all)
  return(results)
}

### unregularized group-level weights
#' @param t_vec vector of task membership
#' @param x list of Kernel matrices
#' @param y_mkl MKL label vector
#' @return M MKL matrix
constructGroupMKL <- function(x){
  
  # data features
  names_tasks <- unique(t_vec)
  
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
  M <- as.data.frame(as.matrix(M))
  
  # reorder matrix to match original indices
  indices_tasks <- as.numeric(unlist(indices_tasks))
  M <- cbind(M, indices_tasks)
  M <- arrange(M, indices_tasks) %>%
    select(-indices_tasks)
  M <- as.matrix(M)
  
  return(M)
}

### WIP: blockwise MTMKL
#' @param t_vec vector of task membership
#' @param x list of Kernel matrices
#' @param y_mkl MKL label vector
#' @return M MKL matrix
constructBlockMKL <- function(x){
  
  # data features
  names_tasks <- unique(t_vec)
  
  # generate a named list of task indices
  indices_tasks <- vector(mode = "list", length = length(names_tasks))
  for (i in 1:length(names_tasks)) {
    indices_tasks[[i]] <- which(t_vec == names_tasks[[i]])
  }
  names(indices_tasks) <- names_tasks
  
  # construct a task-wise block diagonal matrix first
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
  
  # get new task indices
  indices_tasks_new <- list()
  for (i in 1:length(names_tasks)) {
    l <- length(indices_tasks[[i]]) # length of index vector
    k <- length(unlist(indices_tasks[(1:i)-1])) # length of all previous tasks index vectors
    t <- c(k+(1:l)) # task block coordinates
    indices_tasks_new[[i]] <- t
  }
  names(indices_tasks_new) <- names_tasks
  
  list_2tasks <- expand.grid(names_tasks, names_tasks,
                             stringsAsFactors = FALSE) %>%
    as.data.frame() %>%
    split(., seq(nrow(.)), drop = TRUE)
  
  m_2tasks <- list() # list of combined task matrices
  wt <- list() # list of taskwise weight vectors
  for (i in 1:length(list_2tasks)) {
    i1 <- indices_tasks[[list_2tasks[[i]][[1]]]] # indices of first task to grab from
    i2 <- indices_tasks[[list_2tasks[[i]][[2]]]] # indices of second task to grab from
    m <- list(hpo[c(i1,i2), c(i1,i2)], # view matrices for MKL
              M[c(i1,i2), c(i1,i2)]) # be careful if M was stored previously - we want the original instance-level matrix of sequence/structure-based feature RBF kernel fn here
    
    w <- list()
    y_d <- vector()
    
    for (task in 1:length(names_tasks)) {
      y_d <- y_mkl[c(i1,i2)]
      
      # small tasks with similar observations may lead to computationally singular systems
      # if this occurs, return the uniformly weighted kernel matrix
      tryCatch(
        expr = {
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
  }
  
  # name stored weight vector
  names(wt) <- lapply(lapply(list_2tasks, unlist), as.character)
  wt <- as_tibble(wt)
  rownames(wt) <- c("hpo", "structure")
  wt <- t(wt)
  assign("wt", wt, envir = .GlobalEnv)
  
  # create block diagonal matrix
  M <- as.matrix(bdiag(m_tasks))
  
  # populate block diagonal matrix with taskwise weighted MKL matrices
  for (i in 1:length(list_2tasks)){
    j1 <- indices_tasks_new[[list_2tasks[[i]][[1]]]] # indices of first task to store later
    j2 <- indices_tasks_new[[list_2tasks[[i]][[2]]]] # indices of second task to store later
    j <- c(j1,j2)
    
    l1 <- length(j1)
    l2 <- length(j2)
    
    k1 <- 1:l2
    k2 <- (1:l1)+l2
    
    m <- m_2tasks[[i]] 
    
    M[j1,j2] <- m[k1,k2] # TODO more bugfixing
  }
  
  # reorder matrix to match original indices
  M <- as.data.frame(M)
  indices_tasks <- as.numeric(unlist(indices_tasks))
  M <- cbind(M, indices_tasks)
  M <- arrange(M, indices_tasks) %>%
    select(-indices_tasks)
  M <- as.matrix(M)
  
  M <- M+Kt # get product kernel matrix with task similarity for MTMKL
  
  return(M)
} 