# Predict y_hat (outcome) given X (variables), Z (covariates) and object model. [S3 method for class 'bivas']
predict.bivas <- function(object, Z=NULL, X, ...) {

  p <- nrow(object$gamma_pos)
  n <- nrow(X)
  q <- nrow(object$cov)

  if(is.null(Z)) {
    Z <- matrix(1,n,1)
  } else {
    Z <- as.matrix(Z)
    if (nrow(Z) != n)
      stop("Inputs X and Z do not match.")
    Z <- cbind(1,Z)
  }

  if (ncol(Z) != q){
    stop("Inputs arguments object and Z are not compatible")
  }

  cov  <- coef(object,T)$cov
  beta <- coef(object,T)$beta

  y_hat <- Z %*% cov + X %*% beta

  return(y_hat)
}


predict.bivas_mt <- function(object, Z=NULL, X, ...) {
  l  <- length(X)
  K  <- nrow(object$gamma_pos)
  nn <- as.numeric(sapply(X,nrow))

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
  q <- q + 1

  if (nrow(object$cov) != q){
    stop("Inputs arguments object and Z are not compatible")
  }

  cov  <- coef(object,T)$cov
  beta <- coef(object,T)$beta

  y_hat <- vector("list",l)
  for(j in 1:l){
    y_hat[[j]] <- Z[[j]] %*% cov[,j] + X[[j]] %*% beta[,j]
  }

  return(y_hat)
}
