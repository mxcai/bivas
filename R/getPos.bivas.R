# Get variable level posterior of bivas object
getPos <- function(object, ...) UseMethod("getPos")


getPos.bivas <- function(object, weight=T){
  pos <- list()
  if (weight) {
    pos$var_pos <- with(object,(gamma_pos%*%weight) * (eta_pos%*%weight))
    pos$group_pos <- with(object, group_pos%*%weight)
  } else {
    pos$var_pos <- with(object,gamma_pos*eta_pos)
    pos$group_pos <- object$group_pos
  }

  rownames(pos$var_pos) <- object$varname
  rownames(pos$group_pos) <- object$groupname

  return(pos)
}

getPos.bivas_mt <- function(object, weight=T){
  pos <- list()
  if (weight) {
    pos$var_pos <- with(object,arrayXvec(gamma_pos,weight) * drop(eta_pos%*%weight))
    pos$group_pos <- with(object, eta_pos%*%weight)
  } else {
    pos$var_pos <- with(object,arrayXvec_ew(gamma_pos,eta_pos))
    pos$group_pos <- object$eta_pos
  }

  rownames(pos$var_pos) <- object$varname
  rownames(pos$group_pos) <- object$groupname

  return(pos)
}

