# choose best logodds using K-fold cross validation
cv.bivas <- function(y,X,Z=NULL,group,nfolds=10,maxIter=1500,tol=1e-6,sb2,se2,alpha,logodds,mu,alpha_jk,pi_k,
                     verbose=F,coreNum=1){

  n <- nrow(X)
  p <- ncol(X)
  q <- ifelse(is.null(Z),0,ncol(Z))

  group <- as.factor(group)
  glevel <- levels(group)
  K <- length(glevel)

  # initialization of model parameters
  if (missing(logodds)) {
    logodds <- seq(-log10(K),0,length.out = 20)
  }

  ns <- length(logodds)

  if (missing(pi_k)) {
    pi_k <- matrix(runif(K*ns),K,ns)
    pi_k <- pi_k/matrix(colSums(pi_k),K,ns,byrow = TRUE)
  }

  if (missing(alpha)) {
    alpha <- 1/log(p+1)
  }
  alpha <- rep(alpha,ns)

  if (missing(alpha_jk)) {
    alpha_jk <- matrix(runif(p*ns),p,ns)
    alpha_jk <- alpha_jk/matrix(colSums(alpha_jk),p,ns,byrow = TRUE)
  }

  if (missing(mu)) {
    mu <- matrix(rnorm(p*ns),p,ns)
  }

  if (missing(se2)) {
    se2 <- var(y)/2
  }
  se2 <- rep(se2,ns)

  if (missing(sb2)) {
    sb2 <- 1
  }
  sb2 <- rep(sb2,ns)

  #check covariates
  # if(missing(Z)){
  #   Z <- NULL
  # }else{
  #   Z <- cbind(1,Z)
  # }

  # decide the cv assignments
  idx <- ceiling(sample(1:n)/n*nfolds)

  cv.sb2   <- matrix(0,nfolds,ns)
  cv.se2   <- matrix(0,nfolds,ns)
  cv.alpha <- matrix(0,nfolds,ns)
  testErr  <- matrix(0,nfolds,ns)

  # report settings
  message("Info: Number of SNPs: ", p)
  message("Info: Number of groups: ", K)
  message("Info: Number of covariates: ", q)
  message("Info: Sample size: ", n)
  message("Info: Number of candidate hyperparameter settings: ", ns)
  message("Info: Number of cv folds: ", nfolds)

  cat("start cv process......... total",nfolds,"validation sets \n")

  for(i in 1:nfolds) {
    cat(i,"-th validation set... \n")

    Z_train <- Z[idx!=i,]
    X_train <- X[idx!=i,]
    y_train <- y[idx!=i]

    Z_test  <- Z[idx==i,]
    X_test  <- X[idx==i,]
    y_test  <- y[idx==i]


    fit <- bivas(y_train,X_train,Z_train,group,maxIter,tol,sb2,se2,alpha,logodds,mu,alpha_jk,pi_k,verbose,coreNum)

    cv.sb2[i,]   <- fit$sb2
    cv.se2[i,]   <- fit$se2
    cv.alpha[i,] <- fit$alpha
    testErr[i,]  <- colMeans((y_test-predict(fit,Z_test,X_test))^2)

  }

  cvm.sb2   <- colMeans(cv.sb2)
  cvm.se2   <- colMeans(cv.se2)
  cvm.alpha <- colMeans(cv.alpha)

  cvsd.sb2   <- apply(cv.sb2,MARGIN = 2,sd)
  cvsd.se2   <- apply(cv.se2,MARGIN = 2,sd)
  cvsd.alpha <- apply(cv.alpha,MARGIN = 2,sd)

  testErr     <- colMeans(testErr)
  idx.min     <- which.min(testErr)
  logodds.min <- logodds[idx.min]

  cat("cross validation finished, optimal hyperparameter setting: logodds = ", logodds.min,"\n")
  cat("fit the whole data set based on optimal hyperparameter setting \n")

  bivas.fit <- bivas(y,X,Z,group,maxIter,tol,sb2[1],se2[1],alpha[1],logodds=logodds.min,mu[,idx.min],
                     alpha_jk[,idx.min],pi_k[,idx.min],verbose,coreNum)

  cvfit <- list()
  cvfit$logodds    <- logodds
  cvfit$cvm.sb2    <- cvm.sb2
  cvfit$cvm.se2    <- cvm.se2
  cvfit$cvm.alpha  <- cvm.alpha
  cvfit$cvsd.sb2   <- cvsd.sb2
  cvfit$cvsd.se2   <- cvsd.se2
  cvfit$cvsd.alpha <- cvsd.alpha
  cvfit$testErr    <- testErr
  cvfit$bivas.fit  <- bivas.fit

  attr(cvfit,"class") <- "cv.bivas"
  return(cvfit)
}
