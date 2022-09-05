# SCN
# Functional variant prediction for voltage-gated sodium channels
# functions for multiple-kernel learning

### Simple MKL 
#' adapted functions from RMKL package
#' cf. DOI 10.1186/s12859-019-2992-1
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

### SEMKL
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
#' @param hierarchical run hierarchical decomposition, cf. http://doc.ml.tu-berlin.de/mkl_workshop/papers/christian_widmer.pdf
#' @param graph network graph or similarity/distance matrix of dim TxT
#' @return M training MKL matrix
#' cf. https://www.cmpe.boun.edu.tr/~ethem/files/papers/gonen_icml08.pdf
#' cf. DOI 10.1186/1471-2105-11-S8-S5
#' this method decomposes the problem into t*t MKL learning problems, where weights are learned for each combination of tasks (leaves of the hierarchy) which reduces to a quasi-conformal transformation (Amari & Wu, 1998)
constructBlockMKL <- function(matrices, label, tasks, 
                              hierarchical = FALSE, graph){
  # matrices <- M_mkl
  # label <- y_mkl
  # tasks <- t_vec
  # hierarchical <- TRUE
  # graph <- G
  
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
  
  if(hierarchical == FALSE){
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
    
    # preallocate matrix
    M <- matrix(data = 0, nrow = dim, ncol = dim)
    
    # populate matrix with taskwise weighted MKL matrices
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
  }
  
  if(hierarchical == TRUE){
    requireNamespace("dendextend", quietly = TRUE)
    
    if(missing(graph)){
      stop("Provide network graph or similarity/distance matrix of dim TxT")
    }
    
    # convert distance matrix to hierarchy, decompose and save node labels
    graph <- dist(as.matrix(G), diag = F)
    graph <- as.matrix(graph)
    colnames(graph) <- rownames(graph) <- names(G)
    graph <- as.dist(graph)
    graph <- hclust(graph, method = "average")
    graph <- as.dendrogram(graph) 
    tree <- graph # store for later
    graph <- partition_leaves(graph) 
    
    m_decomp <- list()
    wt <- list()
    
    for (i in 1:length(graph)) {
      # subset by indices of node tasks
      ind <- indices_tasks[graph[[i]]] %>% unlist(use.names = F) 
      m <- lapply(matrices, function(x) x[ind, ind])
      y_d <- y[ind]
      
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
      
      # merge node subset matrices by weights and store in order of node subsets
      m_decomp[[i]] <- Reduce(`+`,Map(`*`, w, m))
    }
    
    # compute product kernel matrix from node subset kernel matrices
    M_list <- lapply(1:length(graph), matrix, data = 0, nrow = dim, ncol = dim)
    
    for (i in 1:length(graph)){
      # get node subset kernel matrices and store in list of matrices length of nodes, dim n*n
      jnd <- indices_tasks[graph[[i]]] %>% unlist(use.names = FALSE) 
      M_list[[i]][jnd, jnd] <- M[jnd, jnd] + m_decomp[[i]]
    }
    
    # add task-similarity matrix from earlier
    M_list[[length(M_list)+1]] <- Kt
    
    # use MKL to find best composite kernel matrix and store weights
    delta <- SEMKL.classification(k = M_list,
                                  outcome = y,
                                  penalty = mkl_cost)$gamma
    
    M <- Reduce(`+`,Map(`*`, delta, M_list))
  }

  # store output
  d <- list()
  d$M <- M # training kernel matrix
  d$gamma <- wt # task-wise weights or weights for each node subset if hierarchical decomposition is used
  if(hierarchical == TRUE){
    d$graph <- graph # list of subtrees with labels of leaves for each node
    d$delta <- delta # weight vector for composite kernel matrix
    d$tree <- tree # stored dendrgram object
  } 
  return(d)
} 

#' @param matrices list of Kernel matrices, where the first element is Kt (task-wise similarity)
#' @param label label vector 
#' @param tasks vector t_vec of task membership 
#' @param gamma task-wise weights learned in constructBlockMKL fn
#' @param hierarchical apply weights from hierarchical decomposition?
#' @param graph list of subtree labels from previous function
#' @param delta vector of weights for composite kernel matrix MKL of length nodes+1 (task-wise similarity matrix)
#' @return M MKL for training/testing
#' this function applies weights learned from the training indices and yields a complete MKL matrix for train/test splits
applyBlockMKL <- function(matrices, tasks, gamma,
                          hierarchical = FALSE, graph = NULL, delta = NULL) {
  # debug
  # matrices <- M_mkl
  # tasks <- t_vec
  # gamma <- d$gamma
  # graph <- d$graph
  
  if(missing(graph)){
    stop("Provide list of decomposed node subsets generated by constructBlockMKL fn.")
  }
  
  if(missing(delta)){
    stop("Provide weight vector for composite kernel matrix, should be length nodes+1")
  }
  
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
  
  if (hierarchical == FALSE){
    # list of combined task matrices
    m_2tasks <- list() 
    
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
  }
  
  if(hierarchical == TRUE){
    # list of node subset kernel matrices
    m_decomp <- list() 
    
    for (i in 1:length(graph)) {
      # subset by indices of node tasks
      ind <- indices_tasks[graph[[i]]] %>% unlist(use.names = F) 
      m <- lapply(matrices, function(x) x[ind, ind])
      
      # merge view matrices by weights
      m_decomp[[i]] <- Reduce(`+`,Map(`*`, gamma[[i]], m))
    }
    
    # compute product kernel matrix from node subset kernel matrices
    M_list <- lapply(1:length(graph), matrix, data = 0, nrow = dim, ncol = dim)
    
    for (i in 1:length(graph)){
      # get node subset kernel matrices and store in list of matrices length of nodes, dim n*n
      jnd <- indices_tasks[graph[[i]]] %>% unlist(use.names = FALSE) 
      M_list[[i]][jnd, jnd] <- M[jnd, jnd] + m_decomp[[i]]
    }
    
    # add task-similarity matrix from earlier
    M_list[[length(M_list)+1]] <- Kt
    
    # apply delta weight vector to compute composite kernel matrix
    M <- Reduce(`+`,Map(`*`, delta, M_list))
  }
  return(M)
}