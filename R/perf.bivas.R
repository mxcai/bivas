# Evaluate the performance of bivas model
perf <- function(object, ...) UseMethod("perf")

perf.bivas <- function(y,X,Z=NULL,group,logodds=NULL,nfolds=10,...) {
  n  <- nrow(X)
  p  <- ncol(X)
  q <- ifelse(is.null(Z),0,ncol(Z))
  K  <- length(unique(group))
  ns <- ifelse(is.null(logodds),40,length(logodds))

  # decide the cv assignments
  idx <- ceiling(sample(1:n)/n*nfolds)


  testErr <- rep(0,nfolds)

  # report settings
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


    fit <- bivas(y=y_train,X=X_train,Z=Z_train,group=group,logodds=logodds,verbose=F,...)

    testErr[i]  <- mean((y_test-predict(fit,Z_test,X_test))^2)

  }

  cvm  <- mean(testErr)
  cvsd <- sd(testErr)

  return(list(cvm=cvm,cvsd=cvsd,testErr=testErr))
}



perf.bivas_mt <- function(y,X,Z=NULL,logodds=NULL,nfolds=10,method = c("mse","cor"),...) {
  nn <- sapply(y,length)
  l  <- length(y)
  K  <- ncol(X[[1]])
  p  <- l * K
  q  <- ifelse(is.null(Z),0,ncol(Z))
  ns <- ifelse(is.null(logodds),40,length(logodds))

  # decide the cv assignments
  idx <- list()
  for(j in 1:l){
    idx[[j]] <- ceiling(sample(1:nn[j])/nn[j]*nfolds)
  }



  if("mse" %in% method){
    testErr <- matrix(0,nfolds,l)
  }
  if("auc" %in% method){
    auc <- matrix(0,nfolds,l)
  }
  if("cor" %in% method){
    cor <- matrix(0,nfolds,l)
  }


  # report settings
  message("Info: Number of cv folds: ", nfolds)

  cat("start cv process......... total",nfolds,"validation sets \n")

  for(i in 1:nfolds) {
    cat(i,"-th validation set... \n")

    if(!is.null(Z)){
      Z_train <- mapply(FUN = function(X,idx) X[idx!=i,],Z,idx,SIMPLIFY = F)
    } else {
      Z_train <- NULL
    }
    X_train <- mapply(FUN = function(X,idx) X[idx!=i,],X,idx,SIMPLIFY = F)
    y_train <- mapply(FUN = function(X,idx) X[idx!=i],y,idx,SIMPLIFY = F)

    if(!is.null(Z)){
      Z_test <- mapply(FUN = function(X,idx) X[idx==i,],Z,idx,SIMPLIFY = F)
    } else {
      Z_test <- NULL
    }
    X_test <- mapply(FUN = function(X,idx) X[idx==i,],X,idx,SIMPLIFY = F)
    y_test <- mapply(FUN = function(X,idx) X[idx==i],y,idx,SIMPLIFY = F)

    fit_mt <- bivas_mt(y=y_train,X=X_train,Z=Z_train,logodds=logodds,verbose=F,...)

    y_hat <- predict(fit_mt,Z_test,X_test)

    if("mse" %in% method){
      testErr[i,] <- mapply(FUN = function(y1,y2) mean((y1-as.numeric(y2))^2),y_test,y_hat,SIMPLIFY = T)
    }
    if("auc" %in% method){
      auc[i,] <- mapply(FUN = function(y1,y2) roc(drop(y1),drop(y2))$auc,y_test,y_hat,SIMPLIFY = T)
    }
    if("cor" %in% method){
      cor[i,] <- mapply(FUN = function(y1,y2) cor(y1,as.numeric(y2)),y_test,y_hat,SIMPLIFY = T)
    }

  }

  ret <- list()
  if("mse" %in% method){

    ret[["msem"]] <- colMeans(testErr)
    ret[["msesd"]] <- apply(testErr,2,sd)
    ret[["testErr"]] <- testErr
  }
  if("auc" %in% method){
    ret[["aucm"]] <- colMeans(auc)
    ret[["auc"]] <- auc
  }
  if("cor" %in% method){
    ret[["corm"]] <- colMeans(cor)
    ret[["cor"]] <- cor
  }

  return(ret)

}
