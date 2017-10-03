# Get fdr from bivas object
fdr <- function(object, ...) UseMethod("fdr")

fdr.bivas <- function(object,FDRset=0.1,control="local"){
  p <- nrow(object$mu_pos)
  K <- nrow(object$group_pos)

  ret <- list()
  ret$FDR <- rep(0,p)
  ret$groupFDR <- rep(0,K)

  pos <- getPos(object,T)
  # var_pos <- with(object, eta_pos%*%weight * gamma_pos%*%weight)
  # group_pos <- with(object, group_pos%*%weight)

  if(control=="global"){
    global       <- fdr2FDR(1-pos$var_pos)
    group_global <- fdr2FDR(1-pos$group_pos)

    ret$FDR[which(global<=FDRset)] <- 1
    ret$groupFDR[which(group_global<=FDRset)] <- 1
  }
  if(control=="local"){
    ret$FDR[which((1-pos$var_pos)<=FDRset)] <- 1
    ret$groupFDR[which((1-pos$group_pos)<=FDRset)] <- 1
  }

  names(ret$FDR)   <- object$varname
  names(ret$groupFDR) <- object$groupname

  return(ret)
}

fdr.bivas_mt <- function(object,FDRset=0.1,control="local"){
  K <- nrow(object$mu_pos)
  l <- ncol(object$mu_pos)
  p <- K * l

  ret <- list()
  ret$FDR <- rep(0,p)
  ret$groupFDR <- rep(0,K)

  pos <- getPos(object,T)
  # var_pos <- with(object, eta_pos%*%weight * gamma_pos%*%weight)
  # group_pos <- with(object, group_pos%*%weight)

  if(control=="global"){
    global       <- fdr2FDR(1-pos$var_pos)
    group_global <- fdr2FDR(1-pos$group_pos)

    ret$FDR[which(global<=FDRset)] <- 1
    ret$FDR <- matrix(ret$FDR,K,l)
    ret$groupFDR[which(group_global<=FDRset)] <- 1
  }
  if(control=="local"){
    ret$FDR[which((1-pos$var_pos)<=FDRset)] <- 1
    ret$FDR <- matrix(ret$FDR,K,l)
    ret$groupFDR[which((1-pos$group_pos)<=FDRset)] <- 1
  }

  names(ret$FDR)   <- object$varname
  names(ret$groupFDR) <- object$groupname

  return(ret)
}
