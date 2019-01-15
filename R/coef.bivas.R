# Get coefficients of bivas object
coef.bivas <- function(object, weight=T){
  coef <- list()
  if (weight) {
    coef$cov  <- with(object,cov%*%weight)
    coef$beta <- with(object,(gamma_pos%*%weight) * (eta_pos%*%weight) * (mu_pos%*%weight))
  } else {
    coef$cov  <- object$cov
    coef$beta <- with(object,gamma_pos*eta_pos*mu_pos)
  }

  rownames(coef$beta) <- object$varname
  rownames(coef$cov)  <- object$covname

  return(coef)
}

coef.bivas_mt <- function(object, weight=T){
  coef <- list()
  if(weight){
    coef$cov  <- with(object,matrix(arrayXvec(cov,weight),nrow = nrow(cov)))
    coef$beta <- with(object,arrayXvec(gamma_pos,weight) * drop(eta_pos%*%weight) * arrayXvec(mu_pos,weight))
  } else {
    coef$cov  <- object$cov
    coef$beta <- with(object,arrayXvec_ew(gamma_pos*mu_pos,eta_pos))
  }

  rownames(coef$beta) <- object$varname
  rownames(coef$cov)  <- object$covname

  return(coef)
}
