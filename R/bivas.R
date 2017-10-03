##bi-level VB function with importance sampling on hyperparameters
bivas_omp <- function (y,X,Z=NULL,group,maxIter=1500,tol=1e-6,sb2,se2,alpha,logodds=NULL,mu,alpha_jk,pi_k,verbose=TRUE,coreNum=1) {

  n <- nrow(X)
  p <- ncol(X)

  group <- as.factor(group)
  glevel <- levels(group)
  K <- length(glevel)

  if (!is.null(colnames(X))) {
    varname <- colnames(X[[1]])
  }
  if (!is.null(colnames(Z))) {
    covname <- colnames(Z[[1]])
  }

  # Initialization of model parameters
  if (missing(alpha))   {alpha <- 1/log(p+1)}
  if (missing(se2))     {se2 <- var(y)/2}
  if (missing(sb2))     {sb2 <- se2}


  if (is.null(logodds)) {
    if(length(alpha) == 1 & length(se2) == 1 & length(sb2) == 1) {
      logodds <- seq(-log10(K),log10(K),length.out = 40)
    } else {
      stop("logodds can only be missing when length(alpha) = length(se2) = length(sb2) = 1")
    }
  }

  # Initialize variational parameters
  if (missing(pi_k)) {
    pi_k <- rep.row(1/(1+exp(-logodds)),K)#matrix(runif(K),K,1)
    # pi_k <- pi_k/sum(pi_k)
  }
  pi_k <- as.matrix(pi_k)

  if (missing(alpha_jk)) {
    if(length(alpha)==1){
      alpha_jk <- rep(alpha,p)
    } else{
      alpha_jk <- rep.row(alpha,p)
    }
    # alpha_jk <- rep(alpha,p)#matrix(runif(p),p,1)
    # alpha_jk <- alpha_jk/sum(alpha_jk)
  }
  alpha_jk <- as.matrix(alpha_jk)

  if (missing(mu)) {
    mu <- matrix(0,p,1)#matrix(rnorm(p),p,1)
  }
  mu <- as.matrix(mu)

  # Adjust initial value of alpha
  # alpha   <- alpha * (1+exp(-logodds))



  # Ensure the numbers of candidate hyperparameter settings are consistent
  ns <- max(length(logodds),length(sb2),length(se2),ncol(pi_k),ncol(alpha_jk),ncol(mu))

  if (length(alpha)==1) {alpha <- rep(alpha,ns)}
  if (length(se2)==1)   {se2 <- rep(se2,ns)}
  if (length(sb2)==1)   {sb2 <- rep(sb2,ns)}
  if (ncol(pi_k) == 1) {
    pi_k <- rep.col(pi_k,ns)
  }
  if (ncol(alpha_jk) == 1) {
    alpha_jk <- rep.col(alpha_jk,ns)
  }
  if (ncol(mu) == 1) {
    mu <- rep.col(mu,ns)
  }

  if (length(alpha) != ns | length(se2) != ns | length(sb2) != ns | length(logodds) != ns | ncol(pi_k) != ns | ncol(alpha_jk) != ns | ncol(mu) != ns) {
    stop("initial parameter settings are not consistent")
  }

  if (nrow(pi_k) != K) {
    stop("input pi_k must have row number the same as the number of groups")
  }

  if (nrow(alpha_jk) != p) {
    stop("input alpha_jk must have as many rows as X has columns")
  }

  if (nrow(mu) != p) {
    stop("input mu must have as many rows as X has columns")
  }


  #convert string to numeric type for cpp processing (for convenience of cpp implementation)
  i_group  <- as.numeric(group)
  i_glevel <- c(1:max(i_group))

  #check covariates
  if(is.null(Z)){
    q <- 0
    Z <- matrix(1,n,1)
  }else{
    q <- ncol(Z)
    Z <- cbind(1,Z)
  }

  # reportin settings
  if(verbose){
    message("Info: Number of SNPs: ", p)
    message("Info: Number of groups: ", K)
    message("Info: Number of covariates: ", q)
    message("Info: Sample size: ", n)
    message("Info: Number of candidate hyperparameter settings: ", ns)
  }

  #initialize storage
  logw <- rep(0,ns)
  Lq <- rep(0,ns)

  #############################################################################
  #cpp starts here
  #############################################################################



  result <- outerloop_omp(X, y, Z, i_group, i_glevel, logodds, sb2, se2, alpha, mu, alpha_jk, pi_k, maxIter, tol, verbose, coreNum)



  ################################################################################
  #cpp part stops here
  ################################################################################

  #compute the normalized (approximate) weight.
  w <- with(result,c(normalizelogweights(logw)))

  fit <- list()

  fit$Lq        <- result$logw
  fit$weight    <- w
  fit$gamma_pos <- result$alpha_jk
  fit$eta_pos   <- result$pi_p
  fit$group_pos <- result$pi_k
  fit$mu_pos    <- result$mu

  fit$cov     <- result$cov
  fit$sb2     <- result$sb2
  fit$se2     <- result$se2
  fit$alpha   <- result$alpha
  fit$Lq_List <- result$Lq_List

  fit$logodds <- logodds
  # fit$X       <- X
  # fit$y       <- y
  # fit$Z       <- Z
  fit$group   <- group

  if (!is.null(colnames(X))) {
    fit$varname <- varname
  }
  if (!is.null(colnames(Z))) {
    fit$covname <- covname
  }
  fit$groupname <- glevel

  attr(fit,"class") <- "bivas"
  return(fit)
}

###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

##bi-level VB function with importance sampling on hyperparameters by multi-threading

bivas <- function (y,X,Z=NULL,group,maxIter=1500,tol=1e-6,sb2,se2,alpha,logodds=NULL,mu,alpha_jk,pi_k,verbose=TRUE,coreNum=1) {

  n <- nrow(X)
  p <- ncol(X)

  group <- as.factor(group)
  glevel <- levels(group)
  K <- length(glevel)

  if (!is.null(colnames(X))) {
    varname <- colnames(X)
  }
  if (!is.null(colnames(Z))) {
    covname <- colnames(Z)
  }

  # Initialization of model parameters
  if (missing(alpha))   {alpha <- 1/log(p+1)}
  if (missing(se2))     {se2 <- var(y)/2}
  if (missing(sb2))     {sb2 <- se2}


  if (is.null(logodds)) {
    if(length(alpha) == 1 & length(se2) == 1 & length(sb2) == 1) {
      logodds <- seq(-log10(K),log10(K),length.out = 40)
    } else {
      stop("logodds can only be missing when length(alpha) = length(se2) = length(sb2) = 1")
    }
  }

  # Initialize variational parameters
  if (missing(pi_k)) {
    pi_k <- rep.row(1/(1+exp(-logodds)),K)#matrix(runif(K),K,1)
    # pi_k <- pi_k/sum(pi_k)
  }
  pi_k <- as.matrix(pi_k)

  if (missing(alpha_jk)) {
    if(length(alpha)==1){
      alpha_jk <- rep(alpha,p)
    } else{
      alpha_jk <- rep.row(alpha,p)
    }
    # alpha_jk <- rep(alpha,p)#matrix(runif(p),p,1)
    # alpha_jk <- alpha_jk/sum(alpha_jk)
  }
  alpha_jk <- as.matrix(alpha_jk)

  if (missing(mu)) {
    mu <- matrix(0,p,1)#matrix(rnorm(p),p,1)
  }
  mu <- as.matrix(mu)

  # Adjust initial value of alpha
  # alpha   <- alpha * (1+exp(-logodds))



  # Ensure the numbers of candidate hyperparameter settings are consistent
  ns <- max(length(logodds),length(sb2),length(se2),ncol(pi_k),ncol(alpha_jk),ncol(mu))

  if (length(alpha)==1) {alpha <- rep(alpha,ns)}
  if (length(se2)==1)   {se2 <- rep(se2,ns)}
  if (length(sb2)==1)   {sb2 <- rep(sb2,ns)}
  if (ncol(pi_k) == 1) {
    pi_k <- rep.col(pi_k,ns)
  }
  if (ncol(alpha_jk) == 1) {
    alpha_jk <- rep.col(alpha_jk,ns)
  }
  if (ncol(mu) == 1) {
    mu <- rep.col(mu,ns)
  }

  if (length(alpha) != ns | length(se2) != ns | length(sb2) != ns | length(logodds) != ns | ncol(pi_k) != ns | ncol(alpha_jk) != ns | ncol(mu) != ns) {
    stop("initial parameter settings are not consistent")
  }

  if (nrow(pi_k) != K) {
    stop("input pi_k must have row number the same as the number of groups")
  }

  if (nrow(alpha_jk) != p) {
    stop("input alpha_jk must have as many rows as X has columns")
  }

  if (nrow(mu) != p) {
    stop("input mu must have as many rows as X has columns")
  }


  #convert string to numeric type for cpp processing (for convenience of cpp implementation)
  i_group  <- as.numeric(group)
  i_glevel <- c(1:max(i_group))

  #check covariates
  if(is.null(Z)){
    q <- 0
    Z <- matrix(1,n,1)
  }else{
    q <- ncol(Z)
    Z <- cbind(1,Z)
  }

  # report settings
  message("Info: Number of variables: ", p)
  message("Info: Number of groups: ", K)
  message("Info: Number of covariates: ", q)
  message("Info: Sample size: ", n)
  message("Info: Number of candidate hyperparameter settings: ", ns)
  if(coreNum >= 2){
    verbose = F
  }

  #initialize storage
  logw <- rep(0,ns)
  Lq <- rep(0,ns)

  #############################################################################
  #cpp starts here
  #############################################################################



  result <- outerloop_thread(X, y, Z, i_group, i_glevel, logodds, sb2, se2, alpha, mu, alpha_jk, pi_k, maxIter, tol, verbose, coreNum)



  ################################################################################
  #cpp part stops here
  ################################################################################

  #compute the normalized (approximate) weight.
  w <- with(result,c(normalizelogweights(logw)))

  fit <- list()

  fit$Lq        <- result$logw
  fit$weight    <- w
  fit$gamma_pos <- result$alpha_jk
  fit$eta_pos   <- result$pi_p
  fit$group_pos <- result$pi_k
  fit$mu_pos    <- result$mu

  fit$cov     <- result$cov
  fit$sb2     <- result$sb2
  fit$se2     <- result$se2
  fit$alpha   <- result$alpha
  fit$Lq_List <- result$Lq_List

  fit$logodds <- logodds
  # fit$X       <- X
  # fit$y       <- y
  # fit$Z       <- Z
  fit$group   <- group

  if (!is.null(colnames(X))) {
    fit$varname <- varname
  }
  if (!is.null(colnames(Z))) {
    fit$covname <- covname
  }
  fit$groupname <- glevel

  attr(fit,"class") <- "bivas"
  return(fit)
}



###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

##bi-level VB for multi-task learning

bivas_mt <- function (y,X,Z=NULL,maxIter=1500,tol=1e-6,sb2,se2,alpha,logodds=NULL,mu,alpha_jk,pi_k,verbose=TRUE,coreNum=1) {

  if(!is.list(y)){
    stop("input y must be a list.")
  }
  if(!is.list(X)){
    stop("input X must be a list.")
  }

  nn <- sapply(y,length)

  if(any(nn!=sapply(X,nrow))){
    stop("each input of y must have as many elements as X has rows")
  }

  K <- sapply(X,ncol)

  if(any(K[1]!=K)){
    stop("each input of X must have the same number of columns")
  } else {
    K <- K[1]
  }

  l <- length(y)


  if (!is.null(colnames(X[[1]]))) {
    varname <- colnames(X[[1]])
  }
  if (!is.null(colnames(Z[[1]]))) {
    covname <- colnames(Z[[1]])
  }

  # Initialization of model parameters
  if (missing(alpha))   {alpha <- 1/log(K*l+1)}
  if (missing(se2))     {se2 <- matrix(sapply(y,var)/2,l,1)}
  if (missing(sb2))     {sb2 <- matrix(se2,l,1)}


  if (is.null(logodds)) {
    if(length(alpha) == 1 & length(se2) == l & length(sb2) == l) {
      logodds <- seq(-log10(K),log10(K),length.out = 20)
    } else {
      stop("logodds can only be missing when length(alpha) = 1, length(se2) = length(sb2) = #task")
    }
  }

  # Initialize variational parameters
  if (missing(pi_k)) {
    pi_k <- rep.row(1/(1+exp(-logodds)),K)#matrix(runif(K),K,1)
    # pi_k <- pi_k/sum(pi_k)
  }
  pi_k <- as.matrix(pi_k)

  if (missing(alpha_jk)) {
    if(length(alpha)==1){
      alpha <- alpha / (1/(1+exp(-logodds)))
      alpha[alpha>=1] <- 0.9
      alpha_jk <- array(alpha,dim=c(K,l,1))
    } else{
      alpha_jk <- array(rep(alpha,rep(K*l,length(alpha))),c(K,l,length(alpha)))
    }
  }

  if (missing(mu)) {
    mu <- array(0,dim=c(K,l,1))
  }

  # Adjust initial value of alpha
  # alpha   <- alpha * (1+exp(-logodds))



  # Ensure the numbers of candidate hyperparameter settings are consistent
  ns <- max(length(logodds),ncol(sb2),ncol(se2),ncol(pi_k),dim(alpha_jk)[3],dim(mu)[3])

  if (length(alpha)==1) {alpha <- rep(alpha,ns)}
  if (length(se2)==l)   {se2 <- rep.col(se2,ns)}
  if (length(sb2)==l)   {sb2 <- rep.col(sb2,ns)}
  if (ncol(pi_k) == 1) {
    pi_k <- rep.col(pi_k,ns)
  }
  if (dim(alpha_jk)[3] == 1) {
    alpha_jk <- array(rep(alpha_jk,ns),c(K,l,ns))
  }
  if (dim(mu)[3] == 1) {
    mu <- array(rep(mu,ns),c(K,l,ns))
  }

  if (length(alpha) != ns | ncol(se2) != ns | ncol(sb2) != ns | length(logodds) != ns | ncol(pi_k) != ns | dim(alpha_jk)[3] != ns | dim(mu)[3] != ns) {
    stop("initial parameter settings are not consistent")
  }

  if (nrow(pi_k) != K) {
    stop("input pi_k must have row number the same as the number of groups")
  }

  if (nrow(alpha_jk) != K | ncol(alpha_jk) != l) {
    stop("input alpha_jk has inconsistent dimensionality")
  }

  if (nrow(mu) != K | ncol(mu) != l) {
    stop("input mu has inconsistent dimensionality")
  }


  #check covariates
  if(is.null(Z)){
    q <- 0
    Z <- vector("list",l)
    for(j in 1:l){
      Z[[j]] <- matrix(1,nn[j],1)
    }
  }else{
    q <- sapply(Z,ncol)
    if(any(q[1]!=q)){
      stop("each input of Z must have the same number of columns")
    } else {
      q <- q[1]
    }
    Z <- lapply(Z,function(x) cbind(1,x))
  }

  y <- lapply(y,as.matrix)
  X <- lapply(X,as.matrix)

  # report settings
  message("Info: Number of variables: ", K)
  message("Info: Number of tasks: ", l)
  message("Info: Number of covariates: ", q)
  # message("Info: Sample size: ", nn)
  message("Info: Number of candidate hyperparameter settings: ", ns)
  if(coreNum >= 2){
    verbose = F
  }

  #initialize storage
  logw <- rep(0,ns)
  Lq <- rep(0,ns)

  #############################################################################
  #cpp starts here
  #############################################################################



  result <- outerloop_mt_thread(X, y, Z, logodds, sb2, se2, alpha, mu, alpha_jk, pi_k, nn, K, q, maxIter, tol, verbose, coreNum)



  ################################################################################
  #cpp part stops here
  ################################################################################

  #compute the normalized (approximate) weight.
  w <- with(result,c(normalizelogweights(logw)))

  fit <- list()

  fit$Lq        <- result$logw
  fit$weight    <- w
  fit$gamma_pos <- result$alpha_jk
  fit$eta_pos <- result$pi_k
  fit$mu_pos    <- result$mu

  fit$cov     <- result$cov
  fit$sb2     <- result$sb2
  fit$se2     <- result$se2
  fit$alpha   <- result$alpha
  fit$Lq_List <- result$Lq_List

  fit$logodds <- logodds
  # fit$X       <- X
  # fit$y       <- y
  # fit$Z       <- Z

  if (!is.null(colnames(X[[1]]))) {
    fit$varname <- varname
  }
  if (!is.null(colnames(Z[[1]]))) {
    fit$covname <- covname
  }

  attr(fit,"class") <- "bivas_mt"
  return(fit)
}

