\name{bivas-package}
\alias{bivas-package}
\docType{package}
\title{
bivas (A scalable Bayesian method for bi-level variable selection)
}
\description{
This package provides functions for fitting bivas,
a statistical approach for bi-level variable selection.
}
\details{
\tabular{ll}{
Package: \tab bivas\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2017-09-04\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}

}
\author{
Mingxuan Cai

Maintainer: Mingxuan Cai <12250848@life.hkbu.edu.hk>
}
\references{

}
\keyword{ package }
\seealso{
 \code{\link{coef.bivas}}, \code{\link{predict.bivas}}, \code{\link{perf.bivas}}, \code{\link{getPos.bivas}}, \code{\link{fdr.bivas}}.
}
\examples{
library(bivas)

fit_bivas <- bivas(y=y,X=X,group=g)
}
