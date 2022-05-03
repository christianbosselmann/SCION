# SCN
# Functional variant prediction for voltage-gated sodium channels
# RMTL helper functions based on github.com/transbioZI/RMTL

#' @params mat training kernel matrix
#' @params y label vector of training kernel matrix, binary classification (-1;1)
#' @params t task vector of training kernel matrix
#' @params G network graph or similarity matrix of dim TxT
#' @description wrapper; this function takes a split kernel matrix with its task and label vectors and a pre-defined hierarchical task similarity matrix (Widmer et al. 2010). The task similarity is then refined via RMTL.
#' @return d list of input data and RMTL task-similarity matrix d$T
constructRMTL <- function(mat, y, t, G) {
  # mat = M[i,i]
  # y = y[i]
  # t = t_vec[i]
  # G = G
  
  d <- list()
  
  dim <- dim(mat)[1]
  
  # convert label into -1;1 if needed
  if(!is.numeric(y)) {
    y <- as.numeric(y) 
    y[y == 1] <- 1
    y[y == 2] <- -1
  }
  
  # split label vector into a task-wise list of labels
  y <- mat %>%
    as_tibble() %>%
    cbind(t) %>%
    cbind(y) %>%
    group_by(t) %>%
    group_split()
  y <- lapply(y, function(x) x[,dim+2, drop = FALSE])
  y <- lapply(y, as.vector)
  
  # split kernel matrix into a task-wise list of matrix subsets
  mat <- mat %>%
    as_tibble() %>%
    cbind(t) %>%
    group_by(t) %>%
    group_split()
  mat <- lapply(mat, function(x) x[,1:(dim-2), drop = FALSE])
  mat <- lapply(mat, as.matrix)
  
  # note: previous bug was caused by incorrect row/column order - has to be the same as the order of task-wise lists of matrices and labels
  G <- as.data.frame(G)
  rownames(G) <- colnames(G)
  manual_order <- c("SCN10A", "SCN11A", "SCN1A", "SCN2A", "SCN3A", "SCN4A", "SCN5A", "SCN8A", "SCN9A")
  G <- G[, manual_order]
  G <- G[manual_order, ]
  
  # run RMTL
  d$X <- mat
  d$Y <- y
  d$G <- as.matrix(G)
  
  d$cvfit <- cvMTL(d$X, d$Y, 
                   type = "Classification", Regularization = "Graph", 
                   G = d$G, 
                   Lam1_seq = 10^seq(2,-6, -1),  Lam2 = 2, 
                   opts = list(init = 0, tol = 10^-3, maxIter = 1000), 
                   nfolds = 2, stratify = FALSE, parallel = TRUE, ncores = 4)
  
  d$m <- MTL(d$X, d$Y, 
             type = "Classification", Regularization = "Graph", 
             Lam1 = d$cvfit$Lam1.min, Lam1_seq = d$cvfit$Lam1_seq, 
             G = d$G)
  
  # task correlation matrix stored in d$m$W and named as input graph
  d$T <- cor(d$m$W) 
  d$T <- kernelCentering(d$T)
  d$T <- as.data.frame(d$T)
  rownames(d$T) <- colnames(d$G)
  colnames(d$T) <- colnames(d$G)
  
  return(d)
}

#' @params mat kernel matrix
#' @params t task vector
#' @params d returned from constructRMTL function, incl. RMTL task similarity matrix d$T
#' @description wrapper; this function takes a refined RMTL task similarity matrix, transforms it to a kernel matrix Kt and computes a product kernel with the input kernel matrix.
#' @return M product global kernel matrix, to be split into train and test
applyRMTL <- function(mat, t, d){
  G <- d$T
  
  G <- G %>%
    rownames_to_column() %>%
    pivot_longer(cols = -c(1)) %>%
    dplyr::rename(k = rowname, l = name)
  
  m <- matrix(data = NA, nrow = dim(M)[1], ncol = dim(M)[1])
  colnames(m) <- t
  rownames(m) <- t
  
  for(i in 1:nrow(G)){
    m[rownames(m) == G[[i,2]], colnames(m) == G[[i,1]]] <- G[[i,3]]
  }
  
  M <- M*m
  return(M)
}

# core function from RMTL package
MTL <- function(X, Y, type="Classification", Regularization="Graph",
                Lam1=0.1, Lam1_seq=NULL, Lam2=0,
                opts=list(init=0,  tol=10^-3,
                          maxIter=1000), G=NULL, k=2){
  
  # X=cv_Xtr
  # Y=cv_Ytr
  # type=type
  # Regularization=Regularization
  # Lam1=Lam1_seq[p_idx]
  # Lam2=Lam2
  # opts=opt
  # k=k
  # G=G
  
  #test validity of input data
  if (!missing(X) & !missing(Y)){
    if (all(sapply(X, class)!="matrix")){
      X <- lapply(X, function(x){as.matrix(x)})
    }
    if (all(sapply(Y, class)!="matrix")){
      Y <- lapply(Y, function(x){as.matrix(x)})
    }
  }else{
    stop("data X or Y doesnot exists")
  }
  
  #test the validity of problem type
  if(type=="Classification"){
    method <- "LR"
  }else if(type=="Regression"){
    method <- "LS"
  }else{
    stop("neither Regression or Classification")
  }
  
  #test the validity of regularization 
  allRegularizations <- c("L21", "Lasso", "Graph", "CMTL", "Trace")
  if (is.element(Regularization, allRegularizations)){
    method <- paste0(method, "_", Regularization)
  }else{
    stop("Regularization is not recognizable")}
  
  #test validity of Lam1 and Lam2
  if (Lam1<0) {stop("Lam1 must be positive")}
  if (Lam2<0) {stop("Lam2 must be positive")}
  
  #collect arguments 
  args <- list(X=X, Y=Y, lam1=Lam1, lam2=Lam2, opts=opts)
  
  if (Regularization=="CMTL"){
    if (k>0){
      args$k <- k
    }else(stop("for CMTL, k must be positive interger"))
  }
  if (Regularization=="Graph"){
    if(!is.null(G)){
      args$G <- G
    }else{stop("graph matrix G is not provided")}
  }
  
  #call solver
  if (any(Regularization==c("L21", "Lasso", "Trace"))){
    #sparse routine
    if( !is.null(Lam1_seq) & length(Lam1_seq)>0){
      #with warm start
      opt <- opts
      for (x in Lam1_seq){
        args$lam1 <- x
        m <- do.call(method, args)
        opt$init <- 1
        opt$W0 <- m$W
        opt$C0 <- m$C
        args$opts <- opt
        if (x<=Lam1) break
      }
    } else {
      #without warm start
      m <- do.call(method, args)
    }
  } else if(any(Regularization==c("Graph", "CMTL"))){
    m <- do.call(LR_Graph, args) # TODO hack to make sure furrr finds LR_Graph fn
    # m <- do.call(method, args)
     }
  
  m$call <- match.call()
  m$Lam1 <- args$lam1
  m$Lam2 <- args$Lam2
  m$opts <- args$opts
  m$dim <- sapply(X, function(x)dim(x))
  m$type=type
  m$Regularization=Regularization
  m$method=method
  class(m) <- "MTL"
  return(m)
}

# core function from RMTL package
predict.MTL <- function(object, newX=NULL, ...){
  if(!is.null(newX)){
    task_num <- length(newX)
    score <- lapply(c(1:task_num), function(x)
      newX[[x]] %*% object$W[,x] + object$C[x])
    if (object$type=="Classification"){
      y <- lapply(c(1:task_num),
                  function(x) exp(score[[x]]))
      y <- lapply(y, function(x) x/(1+x))
    }else if (object$type=="Regression"){
      y <- score
    }
    return(y)
  }else{stop("no new data (X) is provided")}
}

# core function from RMTL package
calcError <- function(m, newX=NULL, newY=NULL){
  if(class(m)!="MTL"){
    stop("The first arguement is not a MTL model")}
  if(!is.null(newX) & !is.null(newY)){
    task_num <- length(newY)
    yhat <- predict.MTL(m,newX)
    if(m$type=="Classification"){
      residue <- lapply(1:task_num, function(x)
        newY[[x]]-(round(yhat[[x]])-0.5)*2)
      error <- sapply(residue, function(x){sum(x!=0)/length(x)})
    }else if(m$type=="Regression"){
      error <- sapply(1:task_num, function(x)
        mean((newY[[x]]-yhat[[x]])^2))
    }
    return(mean(error))
  }else{stop(" no new data (X or Y) are provided ")}
}

# core function from RMTL package
cvMTL <- function(X, Y, type="Classification", Regularization="Graph",
                  Lam1_seq=10^seq(0,-6, -1), Lam2=0, G=G, k=2,
                  opts=list(init=0, tol=10^-3, maxIter=1000),
                  stratify=FALSE, nfolds=5, ncores=2, parallel=FALSE){

  # X <- d$X
  # Y <- d$Y
  # G <- d$G
  # type="Classification"
  # Regularization="Graph"
  # Lam1_seq=10^seq(1,-4, -1)
  # Lam2=0
  # k=2
  # opts=list(init=0, tol=10^-3, maxIter=1000)
  # stratify=FALSE
  # nfolds=5
  # ncores=2
  # parallel=FALSE
  
  #test validity of input data
  if (!missing(X) & !missing(Y)){
    if (all(sapply(X, class)!="matrix")){
      X <- lapply(X, function(x){as.matrix(x)})
    }
    if (all(sapply(Y, class)!="matrix")){
      Y <- lapply(Y, function(x){as.matrix(x)})
    }
  }else{
    stop("data X or Y doesnot exists")
  }
  task_num <- length(X)
  if(stratify & type=="Regression"){
    stop("stratified CV is not applicable to regression")}
  cvPar <- getCVPartition(Y, nfolds, stratify)
  
  #cv
  if (!parallel){
    cvm <- rep(0, length(Lam1_seq))
    for (i in 1:nfolds){
      cv_Xtr <- lapply(c(1:task_num),
                       function(x) X[[x]][cvPar[[i]][[1]][[x]], , drop = FALSE]) # another small bugfix to avoid errors with task_num = 1
      cv_Ytr <- lapply(c(1:task_num),
                       function(x) Y[[x]][cvPar[[i]][[1]][[x]]])
      cv_Xte <- lapply(c(1:task_num),
                       function(x) X[[x]][cvPar[[i]][[2]][[x]], , drop = FALSE]) # another small bugfix to avoid errors with task_num = 1
      cv_Yte <- lapply(c(1:task_num),
                       function(x) Y[[x]][cvPar[[i]][[2]][[x]]])
      opt <- opts
      for (p_idx in 1:length(Lam1_seq)){
        m <- MTL(X=cv_Xtr, Y=cv_Ytr, type=type,
                 Regularization=Regularization, Lam1=Lam1_seq[p_idx],
                 Lam2=Lam2, opts=opt, k=k, G=G)
        
        #non sparse model training
        if (!is.element(Regularization, c("Graph", "CMTL"))){
          opt$init <- 1
          opt$W0 <- m$W
          opt$C0 <- m$C
        }
        
        cv_err <- calcError(m, newX=cv_Xte, newY=cv_Yte)
        cvm[p_idx] = cvm[p_idx]+cv_err
      }
    }
    cvm = cvm/nfolds
  } else {
    requireNamespace('doParallel', quietly = TRUE)
    requireNamespace('foreach', quietly = TRUE)
    doParallel::registerDoParallel(ncores)
    cvm <- foreach::foreach(i = 1:nfolds, .combine="cbind") %dopar%{
      cv_Xtr <- lapply(c(1:task_num),
                       function(x) X[[x]][cvPar[[i]][[1]][[x]], , drop = FALSE])
      cv_Ytr <- lapply(c(1:task_num),
                       function(x) Y[[x]][cvPar[[i]][[1]][[x]]])
      cv_Xte <- lapply(c(1:task_num),
                       function(x) X[[x]][cvPar[[i]][[2]][[x]], , drop = FALSE])
      cv_Yte <- lapply(c(1:task_num),
                       function(x) Y[[x]][cvPar[[i]][[2]][[x]]])
      opt <- opts
      cvVec=rep(0, length(Lam1_seq))
      for (p_idx in 1: length(Lam1_seq)){
        m <- MTL(X=cv_Xtr, Y=cv_Ytr, type=type,
                 Regularization=Regularization, Lam1=Lam1_seq[p_idx],
                 Lam2=Lam2, opts=opt, k=k, G=G)
        #non sparse model training
        if (!is.element(Regularization, c("Graph", "CMTL")) ){
          opt$init <- 1
          opt$W0 <- m$W
          opt$C0 <- m$C
        }
        cv_err <- calcError(m, newX=cv_Xte, newY=cv_Yte)
        cvVec[p_idx] <- cv_err
      }
      return(cvVec)
    }
    cvm <- rowMeans(cvm)
  }
  
  best_idx <- which(cvm==min(cvm))[1]
  cv <- list(Lam1_seq=Lam1_seq, Lam1.min=Lam1_seq[best_idx],
             Lam2=Lam2, cvm=cvm)
  class(cv) <- "cvMTL"
  return(cv)
}

# core function from RMTL package
getCVPartition <- function(Y, cv_fold, stratify){
  # Y <- Y
  # cv_fold <- nfolds
  # stratify <- FALSE
  
  task_num = length(Y);
  
  randIdx <- lapply(Y, function(x) sample(1:length(x),
                                          length(x), replace = FALSE))        
  cvPar = {};
  for (cv_idx in 1: cv_fold){
    # buid cross validation data splittings for each task.
    cvTrain = {};
    cvTest = {};
    
    #stratified cross validation
    for (t in 1: task_num){
      task_sample_size <- length(Y[[t]]);
      
      if (stratify){ # own fix to avoid error for tasks with sample size 1
        ct <- which(Y[[t]][randIdx[[t]]]==-1);
        cs <- which(Y[[t]][randIdx[[t]]]==1);
        ct_idx <- seq(cv_idx, length(ct), cv_fold);
        cs_idx <- seq(cv_idx, length(cs), cv_fold);
        if(task_sample_size == 1){
          te_idx <- 1
          tr_idx <- 1
        } 
        if(task_sample_size != 1){
          te_idx <- c(ct[ct_idx], cs[cs_idx]);
          tr_idx <- seq(1,task_sample_size)[
            !is.element(1:task_sample_size, te_idx)];
        }
        
      } else { # own fix to avoid error for tasks with sample size 1
        if(task_sample_size == 1){
          te_idx <- 1
          tr_idx <- 1
        } 
        if(task_sample_size != 1){
          te_idx <- seq(cv_idx, task_sample_size, by=cv_fold)
          tr_idx <- seq(1,task_sample_size)[
            !is.element(1:task_sample_size, te_idx)];
        }
        
      }
      
      cvTrain[[t]] = randIdx[[t]][tr_idx]
      cvTest[[t]] = randIdx[[t]][te_idx]
    }
    
    cvPar[[cv_idx]]=list(cvTrain, cvTest);
  }
  return(cvPar)
}

#################################
#least-square solver for regression
#################################
LS_L21 <- function (X, Y, lam1, lam2, opts){
  #------------------------------------------------------
  # private functions
  gradVal_eval <- function (W, C){
    r <- lapply(c(1:task_num), function(x)
      LS_grad_eval(W[, x], C[x], X[[x]], Y[[x]]))
    grad_W <- sapply(r, function(x)x[[1]]) + 2* lam2 * W
    grad_C <- sapply(r, function(x)x[[2]])
    funcVal = sum(sapply(r, function(x)x[[3]])) + lam2 * norm(W, 'f')^2
    return(list(grad_W, grad_C, funcVal))
  }
  
  funVal_eval <- function (W, C){
    return(sum(sapply(c(1:task_num), function(x)
      LS_funcVal_eval(W[, x], C[x], X[[x]], Y[[x]])))+
        lam2 * norm(W, 'f')^2)
  }
  
  nonsmooth_eval <- function (W, lam1){
    return(sum(sqrt(rowSums(W^2)))*lam1)
  }
  #------------------------------------------------------
  
  task_num <- length (X);
  dimension = dim(X[[1]])[2];
  Obj <- vector(); 
  
  #initialize a starting point
  if(opts$init==0){
    W0 <- matrix(0, nrow=dimension, ncol=task_num);
    C0 <- rep(0, task_num);
  }else if(opts$init==1){
    W0 <- opts$W0
    C0 <- opts$C0
  }    
  
  bFlag <- 0; 
  Wz <- W0;
  Cz <- C0;
  Wz_old <- W0;
  Cz_old <- C0;
  
  t <- 1;
  t_old <- 0;
  iter <- 0;
  gamma <- 1;
  gamma_inc <- 2;
  
  while (iter < opts$maxIter){
    alpha <- (t_old - 1) /t;
    
    Ws <- (1 + alpha) * Wz - alpha * Wz_old;
    Cs <- (1 + alpha) * Cz - alpha * Cz_old;
    
    # compute function value and gradients of the search point
    r <- gradVal_eval(Ws, Cs);
    gWs <- r[[1]]
    gCs <- r[[2]]
    Fs <- r[[3]]
    
    
    # the Armijo Goldstein line search scheme
    while (TRUE){
      Wzp <- L21_projection(Ws - gWs/gamma, lam1 / gamma);
      Czp <- Cs - gCs/gamma;
      Fzp <- funVal_eval  (Wzp, Czp);
      
      delta_Wzp <- Wzp - Ws;
      delta_Czp <- Czp - Cs;
      nrm_delta_Wzp <- norm(delta_Wzp, 'f')^2;
      nrm_delta_Czp <- sum(delta_Czp * delta_Czp);
      r_sum <- (nrm_delta_Wzp+nrm_delta_Czp)/2;
      
      Fzp_gamma = Fs + sum(delta_Wzp* gWs) + 
        sum(delta_Czp * gCs) + gamma * r_sum
      
      if (r_sum <=1e-20){
        bFlag=1; 
        break;
      }
      
      if (Fzp <= Fzp_gamma) break else {gamma = gamma * gamma_inc}
      
    }
    
    Wz_old = Wz;
    Cz_old = Cz;
    Wz = Wzp;
    Cz = Czp;
    
    Obj = c(Obj, Fzp + nonsmooth_eval(Wz, lam1));
    
    
    #test stop condition.
    if (bFlag) break;
    if (iter>=2){
      if (abs( Obj[length(Obj)] - Obj[length(Obj)-1] ) <= opts$tol)
        break;
    }
    
    iter = iter + 1;
    t_old = t;
    t = 0.5 * (1 + (1+ 4 * t^2)^0.5);
    
  }
  
  W = Wzp;
  C = Czp;
  return(list(W=W, C=C, Obj=Obj))
}


#################################
#logistic regression solver for classification
#################################
LR_L21 <- function (X, Y, lam1, lam2, opts){
  #------------------------------------------------------------
  # private functions
  gradVal_eval <- function (W, C){
    r <- lapply(c(1:task_num),
                function(x)LR_grad_eval( W[, x], C[x], X[[x]], Y[[x]]))
    grad_W <- sapply(r, function(x)x[[1]]) + 2* lam2 * W
    grad_C <- sapply(r, function(x)x[[2]])
    funcVal = sum(sapply(r, function(x)x[[3]])) + lam2 * norm(W, 'f')^2
    return(list(grad_W, grad_C, funcVal))
  }    
  
  funVal_eval <- function (W, C){
    return(sum(sapply(c(1:task_num),
                      function(x)LR_funcVal_eval(W[, x], C[x], X[[x]], Y[[x]]))) +
             lam2 * norm(W, 'f')^2)
  }
  
  nonsmooth_eval <- function (W, lam1){
    return(non_smooth_value = sum(sqrt(rowSums(W^2)))*lam1)
  }
  #------------------------------------------------------------
  
  #main algorithm
  task_num <- length (X);
  dimension = dim(X[[1]])[2];
  subjects <- dim(X[[1]])[1];
  Obj <- vector(); 
  
  #initialize a starting point
  if(opts$init==0){
    W0 <- matrix(0, nrow=dimension, ncol=task_num);
    C0 <- rep(0, task_num);
  }else if(opts$init==1){
    W0 <- opts$W0
    C0 <- opts$C0
  }    
  
  bFlag <- 0; 
  Wz <- W0;
  Cz <- C0;
  Wz_old <- W0;
  Cz_old <- C0;
  
  t <- 1;
  t_old <- 0;
  iter <- 0;
  gamma <- 1;
  gamma_inc <- 2;
  
  while (iter < opts$maxIter){
    alpha <- (t_old - 1) /t;
    
    Ws <- (1 + alpha) * Wz - alpha * Wz_old;
    Cs <- (1 + alpha) * Cz - alpha * Cz_old;
    
    # compute function value and gradients of the search point
    r <- gradVal_eval(Ws, Cs);
    gWs <- r[[1]]
    gCs <- r[[2]]
    Fs <- r[[3]]
    
    
    # the Armijo Goldstein line search scheme
    while (TRUE){
      Wzp <- L21_projection(Ws - gWs/gamma, lam1 / gamma);
      Czp <- Cs - gCs/gamma;
      Fzp <- funVal_eval  (Wzp, Czp);
      
      delta_Wzp <- Wzp - Ws;
      delta_Czp <- Czp - Cs;
      nrm_delta_Wzp <- norm(delta_Wzp, 'f')^2;
      nrm_delta_Czp <- sum(delta_Czp * delta_Czp);
      r_sum <- (nrm_delta_Wzp+nrm_delta_Czp)/2;
      
      Fzp_gamma = Fs + sum(delta_Wzp* gWs) + 
        sum(delta_Czp * gCs) + gamma/2 * nrm_delta_Wzp +
        gamma/2 * nrm_delta_Czp;
      
      if (r_sum <=1e-20){
        bFlag=1; 
        break;
      }
      
      if (Fzp <= Fzp_gamma) break else {gamma = gamma * gamma_inc}
      
    }
    
    Wz_old = Wz;
    Cz_old = Cz;
    Wz = Wzp;
    Cz = Czp;
    
    Obj = c(Obj, Fzp + nonsmooth_eval(Wz, lam1));
    
    
    #test stop condition.
    if (bFlag) break;
    if (iter>=2){
      if (abs( Obj[length(Obj)] - Obj[length(Obj)-1] ) <= opts$tol)
        break;
    }
    
    iter = iter + 1;
    t_old = t;
    t = 0.5 * (1 + (1+ 4 * t^2)^0.5);
    
  }
  
  
  return(list(W=Wzp, C=Czp, Obj=Obj))
}

#################################
#least-square solver for regression
#################################
LS_Lasso <- function (X, Y, lam1, lam2, opts){
  #------------------------------------------------------
  # private functions
  gradVal_eval <- function (W, C){
    r <- lapply(c(1:task_num), function(x)
      LS_grad_eval(W[, x], C[x], X[[x]], Y[[x]]))
    grad_W <- sapply(r, function(x)x[[1]]) + 2* lam2 * W
    grad_C <- sapply(r, function(x)x[[2]])
    funcVal = sum(sapply(r, function(x)x[[3]])) + lam2 * norm(W, 'f')^2
    return(list(grad_W, grad_C, funcVal))
  }
  
  funVal_eval <- function (W, C){
    return(sum(sapply(c(1:task_num), function(x)
      LS_funcVal_eval(W[, x], C[x], X[[x]], Y[[x]]))) +
        lam2 * norm(W, 'f')^2)
  }
  
  nonsmooth_eval <- function (W, lam1){
    return(lam1*sum(abs(W)))
  }
  #-------------------------------------------------------    
  
  # Main algorithm
  task_num <- length (X);
  dimension = dim(X[[1]])[2];
  Obj <- vector(); 
  
  #initialize a starting point
  if(opts$init==0){
    W0 <- matrix(0, nrow=dimension, ncol=task_num);
    C0 <- rep(0, task_num);
  }else if(opts$init==1){
    W0 <- opts$W0
    C0 <- opts$C0
  }    
  
  bFlag <- 0; 
  Wz <- W0;
  Cz <- C0;
  Wz_old <- W0;
  Cz_old <- C0;
  
  t <- 1;
  t_old <- 0;
  iter <- 0;
  gamma <- 1;
  gamma_inc <- 2;
  
  while (iter < opts$maxIter){
    alpha <- (t_old - 1) /t;
    
    Ws <- (1 + alpha) * Wz - alpha * Wz_old;
    Cs <- (1 + alpha) * Cz - alpha * Cz_old;
    
    # compute function value and gradients of the search point
    r <- gradVal_eval(Ws, Cs);
    gWs <- r[[1]]
    gCs <- r[[2]]
    Fs <- r[[3]]
    
    
    # the Armijo Goldstein line search scheme
    while (TRUE){
      Wzp <- l1_projection(Ws - gWs/gamma, lam1 / gamma);
      Czp <- Cs - gCs/gamma;
      Fzp <- funVal_eval(Wzp, Czp);
      
      delta_Wzp <- Wzp - Ws;
      delta_Czp <- Czp - Cs;
      nrm_delta_Wzp <- norm(delta_Wzp, 'f')^2;
      nrm_delta_Czp <- sum(delta_Czp * delta_Czp);
      r_sum <- (nrm_delta_Wzp+nrm_delta_Czp)/2;
      
      Fzp_gamma = Fs + sum(delta_Wzp* gWs) + 
        sum(delta_Czp * gCs) + gamma * r_sum
      
      if (r_sum <=1e-20){
        bFlag=1; 
        break;
      }
      
      if (Fzp <= Fzp_gamma) break else {gamma = gamma * gamma_inc}
    }
    
    Wz_old = Wz;
    Cz_old = Cz;
    Wz = Wzp;
    Cz = Czp;
    Obj = c(Obj, Fzp + nonsmooth_eval(Wz, lam1));
    
    
    #test stop condition.
    if (bFlag) break;
    if (iter>=2){
      if (abs( Obj[length(Obj)] - Obj[length(Obj)-1] ) <= opts$tol)
        break;
    }
    
    iter = iter + 1;
    t_old = t;
    t = 0.5 * (1 + (1+ 4 * t^2)^0.5);
    
  }
  
  W = Wzp;
  C = Czp;
  return(list(W=W, C=C, Obj=Obj))
}




LR_Lasso <- function (X, Y, lam1, lam2, opts){
  #------------------------------------------------------------
  # private functions
  
  gradVal_eval <- function (W, C){
    r <- lapply(c(1:task_num),
                function(x)LR_grad_eval( W[, x], C[x], X[[x]], Y[[x]]))
    grad_W <- sapply(r, function(x)x[[1]]) + 2* lam2 * W
    grad_C <- sapply(r, function(x)x[[2]])
    funcVal = sum(sapply(r, function(x)x[[3]])) + lam2 * norm(W, 'f')^2
    return(list(grad_W, grad_C, funcVal))
  }    
  
  funVal_eval <- function (W, C){
    return(sum(sapply(c(1:task_num),
                      function(x)LR_funcVal_eval(W[, x], C[x], X[[x]], Y[[x]]))) +
             lam2 * norm(W, 'f')^2)
  }
  
  nonsmooth_eval <- function (W, lam1){
    return(lam1*sum(abs(W)))
  }
  
  #------------------------------------------------------------
  
  task_num <- length (X);
  dimension = dim(X[[1]])[2];
  subjects <- dim(X[[1]])[1];
  Obj <- vector(); 
  
  #initialize a starting point
  if(opts$init==0){
    W0 <- matrix(0, nrow=dimension, ncol=task_num);
    C0 <- rep(0, task_num);
  }else if(opts$init==1){
    W0 <- opts$W0
    C0 <- opts$C0
  }
  
  bFlag <- 0; 
  Wz <- W0;
  Cz <- C0;
  Wz_old <- W0;
  Cz_old <- C0;
  
  t <- 1;
  t_old <- 0;
  iter <- 0;
  gamma <- 1;
  gamma_inc <- 2;
  
  while (iter < opts$maxIter){
    alpha <- (t_old - 1) /t;
    
    Ws <- (1 + alpha) * Wz - alpha * Wz_old;
    Cs <- (1 + alpha) * Cz - alpha * Cz_old;
    
    # compute function value and gradients of the search point
    r <- gradVal_eval(Ws, Cs);
    gWs <- r[[1]]
    gCs <- r[[2]]
    Fs <- r[[3]]
    
    
    # the Armijo Goldstein line search scheme
    while (TRUE){
      Wzp <- l1_projection(Ws - gWs/gamma, lam1 / gamma);
      Czp <- Cs - gCs/gamma;
      Fzp <- funVal_eval  (Wzp, Czp);
      
      delta_Wzp <- Wzp - Ws;
      delta_Czp <- Czp - Cs;
      nrm_delta_Wzp <- norm(delta_Wzp, 'f')^2;
      nrm_delta_Czp <- sum(delta_Czp * delta_Czp);
      r_sum <- (nrm_delta_Wzp+nrm_delta_Czp)/2;
      
      Fzp_gamma = Fs + sum(delta_Wzp* gWs) + 
        sum(delta_Czp * gCs) + gamma/2 * nrm_delta_Wzp +
        gamma/2 * nrm_delta_Czp;
      
      if (r_sum <=1e-20){
        bFlag=1; 
        break;
      }
      if (Fzp <= Fzp_gamma) break else {gamma = gamma * gamma_inc}
      
    }
    
    Wz_old = Wz;
    Cz_old = Cz;
    Wz = Wzp;
    Cz = Czp;
    Obj = c(Obj, Fzp + nonsmooth_eval(Wz, lam1));
    
    #test stop condition.
    if (bFlag) break;
    if (iter>=2){
      if (abs( Obj[length(Obj)] - Obj[length(Obj)-1] ) <= opts$tol)
        break;
    }
    iter = iter + 1;
    t_old = t;
    t = 0.5 * (1 + (1+ 4 * t^2)^0.5);
  }
  
  W = Wzp;
  C = Czp;
  
  return(list(W=W, C=C, Obj=Obj))
}

#################################
#least-square solver for regression
#################################
LS_Graph <- function (X, Y, G, lam1, lam2, opts){
  #------------------------------------------------------
  # private functions
  gradVal_eval <- function (W, C){
    r <- lapply(c(1:task_num), function(x)
      LS_grad_eval(W[, x], C[x], X[[x]], Y[[x]]))
    grad_W <- sapply(r, function(x)x[[1]]) + 2*lam1*W %*% GGt + 2* lam2*W
    grad_C <- sapply(r, function(x)x[[2]])
    funcVal <- sum(sapply(r, function(x)x[[3]])) +
      lam1*norm(W%*%G,'f')^2 + lam2*norm(W,'f')^2
    
    return(list(grad_W, grad_C, funcVal))
  }
  
  funVal_eval <- function (W, C){
    return(sum(sapply(c(1:task_num),
                      function(x) LS_funcVal_eval(W[, x], C[x], X[[x]], Y[[x]])))+
             lam1*norm(W%*%G,'f')^2 + lam2*norm(W,'f')^2)
  }
  #-------------------------------------------------------    
  
  # Main algorithm
  task_num <- length (X);
  dimension = dim(X[[1]])[2];
  Obj <- vector(); 
  
  #precomputation
  GGt <- G %*% t(G)
  
  #initialize a starting point
  if(opts$init==0){
    W0 <- matrix(0, nrow=dimension, ncol=task_num);
    C0 <- rep(0, task_num);
  }else if(opts$init==1){
    W0 <- opts$W0
    C0 <- opts$C0
  }    
  
  bFlag <- 0; 
  Wz <- W0;
  Cz <- C0;
  Wz_old <- W0;
  Cz_old <- C0;
  
  t <- 1;
  t_old <- 0;
  iter <- 0;
  gamma <- 1;
  gamma_inc <- 2;
  
  while (iter < opts$maxIter){
    alpha <- (t_old - 1) /t;
    
    Ws <- (1 + alpha) * Wz - alpha * Wz_old;
    Cs <- (1 + alpha) * Cz - alpha * Cz_old;
    
    # compute function value and gradients of the search point
    r <- gradVal_eval(Ws, Cs);
    gWs <- r[[1]]
    gCs <- r[[2]]
    Fs <- r[[3]]
    
    
    # the Armijo Goldstein line search scheme
    while (TRUE){
      Wzp <- Ws - gWs/gamma;
      Czp <- Cs - gCs/gamma;
      Fzp <- funVal_eval  (Wzp, Czp);
      
      delta_Wzp <- Wzp - Ws;
      delta_Czp <- Czp - Cs;
      nrm_delta_Wzp <- norm(delta_Wzp, 'f')^2;
      nrm_delta_Czp <- sum(delta_Czp * delta_Czp);
      r_sum <- (nrm_delta_Wzp+nrm_delta_Czp)/2;
      
      
      Fzp_gamma = Fs + sum(delta_Wzp* gWs) + 
        sum(delta_Czp * gCs) + gamma * r_sum
      
      if (r_sum <=1e-20){
        bFlag=1; 
        break;
      }
      
      if (Fzp <= Fzp_gamma) break else {gamma = gamma * gamma_inc}
      
    }
    
    Wz_old = Wz; 
    Cz_old = Cz;
    Wz = Wzp;
    Cz = Czp;
    
    Obj = c(Obj, Fzp );
    
    #test stop condition.
    if (bFlag) break;
    if (iter>=2){
      if (abs( Obj[length(Obj)] - Obj[length(Obj)-1] ) <= opts$tol)
        break;
    }
    
    iter = iter + 1;
    t_old = t;
    t = 0.5 * (1 + (1+ 4 * t^2)^0.5);
    
  }
  
  W = Wzp;
  C = Czp;
  return(list(W=W, C=C, Obj=Obj))
}

#logistic regression solver for classification
LR_Graph <- function (X, Y, G, lam1, lam2, opts){
  
  gradVal_eval <- function (W, C){
    r <- lapply(c(1:task_num),
                function(x) LR_grad_eval( W[, x], C[x], X[[x]], Y[[x]]))
    grad_W <- sapply(r, function(x)x[[1]])+2*lam1*W %*% GGt + 2* lam2 * W
    grad_C <- sapply(r, function(x)x[[2]])
    funcVal = sum(sapply(r, function(x)x[[3]])) + lam1*norm(W%*%G,'f')^2 +
      lam2 * norm(W,'f')^2
    return(list(grad_W, grad_C, funcVal))
  }    
  
  funVal_eval <- function (W, C){
    return(sum(sapply(c(1:task_num),
                      function(x) LR_funcVal_eval(W[, x], C[x], X[[x]], Y[[x]])))+
             lam1*norm(W%*%G,'f')^2) + lam2 * norm(W,'f')^2
  }
  
  #main algorithm    
  task_num <- length (X);
  dimension = dim(X[[1]])[2];
  subjects <- dim(X[[1]])[1];
  Obj <- vector(); 
  
  #initialize a starting point
  if(opts$init==0){
    W0 <- matrix(0, nrow=dimension, ncol=task_num);
    C0 <- rep(0, task_num);
  }else if(opts$init==1){
    W0 <- opts$W0
    C0 <- opts$C0
  }    
  
  #precomputation
  GGt <- G %*% t(G)
  
  bFlag <- 0; 
  Wz <- W0;
  Cz <- C0;
  Wz_old <- W0;
  Cz_old <- C0;
  
  t <- 1;
  t_old <- 0;
  iter <- 0;
  gamma <- 1;
  gamma_inc <- 2;
  
  while (iter < opts$maxIter){
    alpha <- (t_old - 1) /t;
    
    Ws <- (1 + alpha) * Wz - alpha * Wz_old;
    Cs <- (1 + alpha) * Cz - alpha * Cz_old;
    
    # compute function value and gradients of the search point
    r <- gradVal_eval(Ws, Cs);
    gWs <- r[[1]]
    gCs <- r[[2]]
    Fs <- r[[3]]
    
    # the Armijo Goldstein line search scheme
    while (TRUE){
      Wzp <- Ws - gWs/gamma;
      Czp <- Cs - gCs/gamma;
      Fzp <- funVal_eval  (Wzp, Czp);
      
      delta_Wzp <- Wzp - Ws;
      delta_Czp <- Czp - Cs;
      nrm_delta_Wzp <- norm(delta_Wzp, 'f')^2;
      nrm_delta_Czp <- sum(delta_Czp * delta_Czp);
      r_sum <- (nrm_delta_Wzp+nrm_delta_Czp)/2;
      
      Fzp_gamma = Fs + sum(delta_Wzp* gWs) + 
        sum(delta_Czp * gCs) + gamma * r_sum
      
      if (r_sum <=1e-20){
        bFlag=1; 
        break;
      }
      
      if (Fzp <= Fzp_gamma) break else {gamma = gamma * gamma_inc}
      
    }
    
    Wz_old = Wz;
    Cz_old = Cz;
    Wz = Wzp;
    Cz = Czp;
    
    Obj = c(Obj, Fzp)
    
    #test stop condition.
    if (bFlag) break;
    if (iter>=2){
      if (abs( Obj[length(Obj)] - Obj[length(Obj)-1] ) <= opts$tol)
        break;
    }
    
    iter = iter + 1;
    t_old = t;
    t = 0.5 * (1 + (1+ 4 * t^2)^0.5);
    
  }
  
  W = Wzp;
  C = Czp;
  return(list(W=W, C=C, Obj=Obj))
}

#least-square solver for regression
LS_Trace <- function (X, Y, lam1, lam2, opts){
  
  #------------------------------------------------------
  # private functions
  gradVal_eval <- function (W, C){
    r <- lapply(c(1:task_num), function(x)
      LS_grad_eval(W[, x], C[x], X[[x]], Y[[x]]))
    grad_W <- sapply(r, function(x)x[[1]]) + 2* lam2 * W
    grad_C <- sapply(r, function(x)x[[2]])
    funcVal = sum(sapply(r, function(x)x[[3]])) + lam2 * norm(W, 'f')^2
    return(list(grad_W, grad_C, funcVal))
  }
  
  funVal_eval <- function (W, C){
    return(sum(sapply(c(1:task_num), function(x)
      LS_funcVal_eval(W[, x], C[x], X[[x]], Y[[x]])
      +  lam2 * norm(W, 'f')^2)))
  }
  
  
  #-------------------------------------------------------    
  
  # Main algorithm
  task_num <- length (X);
  dimension = dim(X[[1]])[2];
  Obj <- vector(); 
  
  #initialize a starting point
  if(opts$init==0){
    W0 <- matrix(0, nrow=dimension, ncol=task_num);
    C0 <- rep(0, task_num);
  }else if(opts$init==1){
    W0 <- opts$W0
    C0 <- opts$C0
  }    
  
  bFlag <- 0; 
  Wz <- W0;
  Cz <- C0;
  Wz_old <- W0;
  Cz_old <- C0;
  
  t <- 1;
  t_old <- 0;
  iter <- 0;
  gamma <- 1;
  gamma_inc <- 2;
  
  while (iter < opts$maxIter){
    alpha <- (t_old - 1) /t;
    
    Ws <- (1 + alpha) * Wz - alpha * Wz_old;
    Cs <- (1 + alpha) * Cz - alpha * Cz_old;
    
    # compute function value and gradients of the search point
    r <- gradVal_eval(Ws, Cs);
    gWs <- r[[1]]
    gCs <- r[[2]]
    Fs <- r[[3]]
    
    # the Armijo Goldstein line search scheme
    while (TRUE){
      Wzp <- trace_projection(Ws - gWs/gamma, 2 * lam1 / gamma);
      Wzp_tn <- Wzp[[2]]
      Wzp <- Wzp[[1]]
      Czp <- Cs - gCs/gamma;
      Fzp <- funVal_eval  (Wzp, Czp);
      
      delta_Wzp <- Wzp - Ws;
      delta_Czp <- Czp - Cs;
      nrm_delta_Wzp <- norm(delta_Wzp, 'f')^2;
      nrm_delta_Czp <- sum(delta_Czp * delta_Czp);
      r_sum <- (nrm_delta_Wzp+nrm_delta_Czp)/2;
      
      Fzp_gamma = Fs + sum(delta_Wzp* gWs) + 
        sum(delta_Czp * gCs) + gamma/2 * nrm_delta_Wzp +
        gamma/2 * nrm_delta_Czp;
      
      if (r_sum <=1e-20){
        bFlag=1; 
        break;
      }
      
      if (Fzp <= Fzp_gamma) break else {gamma = gamma * gamma_inc}
      
    }
    
    Wz_old = Wz;
    Cz_old = Cz;
    Wz = Wzp;
    Cz = Czp;
    
    Obj <- c(Obj, Fzp + lam1 * Wzp_tn)
    
    
    #test stop condition.
    if (bFlag) break;
    if (iter>=2){
      if (abs( Obj[length(Obj)] - Obj[length(Obj)-1] ) <= opts$tol)
        break;
    }
    
    iter = iter + 1;
    t_old = t;
    t = 0.5 * (1 + (1+ 4 * t^2)^0.5);
    
  }
  
  W = Wzp;
  C = Czp;
  return(list(W=W, C=C, Obj=Obj))
}

#logistic regression solver for classification
LR_Trace <- function (X, Y, lam1, lam2, opts){
  #------------------------------------------------------------
  # private functions
  gradVal_eval <- function (W, C){
    r <- lapply(c(1:task_num),
                function(x)LR_grad_eval( W[, x], C[x], X[[x]], Y[[x]]))
    grad_W <- sapply(r, function(x)x[[1]])+ 2* lam2 * W
    grad_C <- sapply(r, function(x)x[[2]])
    funcVal = sum(sapply(r, function(x)x[[3]])) + lam2 * norm(W, 'f')^2
    return(list(grad_W, grad_C, funcVal))
  }    
  
  funVal_eval <- function (W, C){
    return(sum(sapply(c(1:task_num),
                      function(x)LR_funcVal_eval(W[, x], C[x], X[[x]], Y[[x]])))+
             lam2 * norm(W, 'f')^2)
  }
  
  #------------------------------------------------------------
  
  task_num <- length (X);
  dimension = dim(X[[1]])[2];
  subjects <- dim(X[[1]])[1];
  Obj <- vector(); 
  
  #initialize a starting point
  if(opts$init==0){
    W0 <- matrix(0, nrow=dimension, ncol=task_num);
    C0 <- rep(0, task_num);
  }else if(opts$init==1){
    W0 <- opts$W0
    C0 <- opts$C0
  }    
  
  bFlag <- 0; 
  Wz <- W0;
  Cz <- C0;
  Wz_old <- W0;
  Cz_old <- C0;
  
  t <- 1;
  t_old <- 0;
  iter <- 0;
  gamma <- 1;
  gamma_inc <- 2;
  
  while (iter < opts$maxIter){
    alpha <- (t_old - 1) /t;
    
    Ws <- (1 + alpha) * Wz - alpha * Wz_old;
    Cs <- (1 + alpha) * Cz - alpha * Cz_old;
    
    # compute function value and gradients of the search point
    r <- gradVal_eval(Ws, Cs);
    gWs <- r[[1]]
    gCs <- r[[2]]
    Fs <- r[[3]]
    
    
    # the Armijo Goldstein line search scheme
    while (TRUE){
      Wzp <- trace_projection(Ws - gWs/gamma, 2 * lam1 / gamma);
      Wzp_tn <- Wzp[[2]]
      Wzp <- Wzp[[1]]
      Czp <- Cs - gCs/gamma;
      Fzp <- funVal_eval  (Wzp, Czp);
      
      delta_Wzp <- Wzp - Ws;
      delta_Czp <- Czp - Cs;
      nrm_delta_Wzp <- norm(delta_Wzp, 'f')^2;
      nrm_delta_Czp <- sum(delta_Czp * delta_Czp);
      r_sum <- (nrm_delta_Wzp+nrm_delta_Czp)/2;
      
      Fzp_gamma = Fs + sum(delta_Wzp* gWs) + 
        +sum(delta_Czp * gCs) + gamma/2 * nrm_delta_Wzp +
        +gamma/2 * nrm_delta_Czp;
      
      if (r_sum <=1e-20){
        bFlag=1; 
        break;
      }
      
      if (Fzp <= Fzp_gamma) break else {gamma = gamma * gamma_inc}
      
    }
    
    Wz_old = Wz;
    Cz_old = Cz;
    Wz = Wzp;
    Cz = Czp;
    
    Obj = c(Obj, Fzp + lam1 * Wzp_tn);
    
    
    #test stop condition.
    if (bFlag) break;
    if (iter>=2){
      if (abs( Obj[length(Obj)] - Obj[length(Obj)-1]) <= opts$tol)
        break;
    }
    
    iter = iter + 1;
    t_old = t;
    t = 0.5 * (1 + (1+ 4 * t^2)^0.5);
    
  }
  
  W = Wzp;
  C = Czp;
  
  return(list(W=W, C=C, Obj=Obj))
}

#gradients and functions evaluation
#logistic regression
LR_grad_eval <- function( w, c, x, y){
  weight <- 1/length(y)
  l <- -y*(x %*% w + c)
  lp <- l
  lp[lp<0] <- 0
  funcVal <- sum(weight * ( log( exp(-lp) +  exp(l-lp) ) + lp ))
  b <- (-weight*y)*(1 - 1/ (1+exp(l)))
  grad_c <- sum(b)
  grad_w <- t(x) %*% b
  return(list(grad_w, grad_c, funcVal))
}
LR_funcVal_eval <- function ( w, c, x, y){
  weight <- 1/length(y)
  l <- -y*(x %*% w + c)
  lp <- l
  lp[lp<0] <- 0
  return(sum(weight * ( log( exp(-lp) +  exp(l-lp) ) + lp )))
}
#Least square
LS_grad_eval <- function( w, c, x, y, xy, xx){
  grad_w <-  t(x) %*% (x %*% w + c - y) / nrow(x)
  grad_c <-  mean(x %*% w + c -y)
  funcVal <- 0.5 * mean((y - x %*% w -c)^2)
  return(list(grad_w, grad_c, funcVal))
}

LS_funcVal_eval <- function ( w, c, x, y){
  return(0.5 * mean((y - x %*% w -c)^2))
}

#Different Projections
l1_projection <- function (W, lambda ){
  p <- abs(W) - lambda
  p[p<0] <- 0
  Wp <- sign(W) * p
  return(Wp)
}

L21_projection <- function (W, lambda ){
  thresfold <- sqrt(rowSums(W^2))
  zeros <- which(thresfold==0)              
  temp <- 1 - lambda/thresfold
  temp <- ifelse(temp<0, 0, temp)
  Wp = matrix(rep(temp, ncol(W)), nrow=length(temp))*W
  Wp[zeros,] <- 0
  return(Wp)
}

trace_projection <- function (W, lambda ){
  requireNamespace('corpcor')
  eigen <- corpcor::fast.svd(W)
  d <- eigen$d
  thresholded_value = d - lambda / 2;
  dp <- thresholded_value * ( thresholded_value > 0 )
  if (length(dp)>1){
    Wp <- eigen$u %*% diag(dp) %*% t(eigen$v)
  }else{
    Wp <- (eigen$u * dp) %*% t(eigen$v)
  }
  return(list(Wp, sum(dp)))
}

singular_projection <- function(Msp, k){
  requireNamespace('MASS')
  eig <- eigen(Msp, symmetric=TRUE)
  EVector <- eig$vector
  EValue  <- eig$values
  DiagSigz <- bsa_ihb(EValue, rep(1,length(EValue)),
                      k, rep(1,length(EValue)))
  DiagSigz <- DiagSigz[[1]]
  Mzp <- EVector %*% diag(DiagSigz) %*% t(EVector)
  Mzp_Pz <- EVector
  Mzp_DiagSigz <- DiagSigz
  return(list(Mzp, Mzp_Pz, Mzp_DiagSigz))
}

bsa_ihb <- function(a,b,r,u){
  # initilization
  break_flag <- 0;
  t_l <- a/b; t_u <- (a - u)/b;
  T <- c(t_l, t_u)
  t_L <- -Inf; t_U <- Inf;
  g_tL <- 0; g_tU <- 0;
  
  iter = 0
  while (length(T)!=0){
    iter <- iter + 1
    g_t <- 0
    t_hat <- stats::median(T)
    
    U <- t_hat < t_u
    M <- (t_u <= t_hat) & (t_hat <= t_l)
    
    if (sum(U)){
      g_t <- g_t + t(b[U]) %*% u[U] 
    }
    
    if (sum(M)){
      g_t <- g_t + sum(b[M] %*% (a[M] - t_hat * b[M]))
    }
    
    if (g_t > r){
      t_L <- t_hat
      T <- T[T > t_hat]
      g_tL = g_t
    } else if (g_t < r){
      t_U <- t_hat
      T <- T[T < t_hat]
      g_tU <- g_t
    } else{
      t_star <- t_hat
      break_flag <- 1
      break
    }
  }
  if (!break_flag){
    t_star <- t_L - (g_tL -r) * (t_U - t_L)/(g_tU - g_tL)     
  }
  temp <- sapply(a - rep(t_star, length(b))*b, function(x)max(0,x))
  x_star <- sapply(temp, function(x) min(x,u))
  return(list(x_star,t_star,iter))
}
