\name{predict.bivas}
\alias{predict.bivas}
\title{Makes prediction based on fitted \code{bivas} object and new data \code{X}.}
\description{Makes prediction based on fitted \code{bivas} object and new data \code{X}.}
\usage{
predict.bivas(object, Z=NULL, X)
}
\arguments{
\item{object}{\code{bivas} object}

\item{X}{a new random effect design matrix of n by p.}

\item{Z}{a new covariate matrix of n by q. Optional. By default, a column of 1's are added.}

}
\details{
  \code{perf.bivas} Makes prediction based on fitted \code{bivas} object and new data \code{X}.
}
\value{
The function returns a vector of predicted values :
  \item{y_hat}{
     The predicted values of y.
}

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Mingxuan Cai
}
\seealso{
 \code{\link{coef.bivas}}, \code{\link{predict.bivas}}, \code{\link{perf.bivas}}, \code{\link{getPos.bivas}}, \code{\link{fdr.bivas}}
}
\examples{
predict.bivas(object, Z=NULL, X)

}
