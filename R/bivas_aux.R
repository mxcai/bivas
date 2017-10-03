

################################################################################################
################################################################################################
################################################################################################
################################################################################################
# Other functions


# ----------------------------------------------------------------------
# Replicate vector x to create an m x n matrix, where m = length(x).
rep.col <- function (x, n)
  matrix(x,length(x),n,byrow = FALSE)

# ----------------------------------------------------------------------
# Replicate vector x to create an m x n matrix, where m = length(x).
rep.row <- function (x, m)
  matrix(x,m,length(x),byrow = TRUE)

# ----------------------------------------------------------------------
# logpexp(x) returns log(1 + exp(x)). The computation is performed in a
# numerically stable manner. For large entries of x, log(1 + exp(x)) is
# effectively the same as x.
logpexp <- function (x) {
  y    <- x
  i    <- which(x < 16)
  y[i] <- log(1 + exp(x[i]))
  return(y)
}

# ----------------------------------------------------------------------
# Use this instead of log(sigmoid(x)) to avoid loss of numerical precision.
logsigmoid <- function (x)
  -logpexp(-x)

# ----------------------------------------------------------------------
# normalizelogweights takes as input an array of unnormalized
# log-probabilities logw and returns normalized probabilities such
# that the sum is equal to 1.
normalizelogweights <- function (logw) {

  # Guard against underflow or overflow by adjusting the
  # log-probabilities so that the largest probability is 1.
  c <- max(logw)
  w <- exp(logw - c)

  # Normalize the probabilities.
  return(w/sum(w))
}

# ----------------------------------------------------------------------
# convert local fdr to global FDR
fdr2FDR <- function(fdr) {
  idx          <- order(fdr)
  sort_fdr     <- fdr[idx]
  FDR          <- cumsum(sort_fdr) / c(1:length(fdr))
  FDR[FDR > 1] <- 1
  FDR[idx]     <- FDR

  return(FDR)
}


# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# arrayXvec
arrayXvec <- function(X,logw) {
  Xw <- apply(X,MARGIN = 2,FUN = function(X,w) X%*%w,w=logw)
}

# ----------------------------------------------------------------------
# arrayXvec_element-wise
arrayXvec_ew <- function(X,y) {
  Xy <- apply(X,MARGIN = 2,FUN = function(X,y) X*y,y=y)
  Xy <- array(t(Xy),dim(X))
}
