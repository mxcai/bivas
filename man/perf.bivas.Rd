\name{perf.bivas}
\alias{perf.bivas}
\title{Evaluate the performance of bivas for a given data using cross-validation.}
\description{Evaluate the performance of bivas for a given data using cross-validation.}
\usage{
perf.bivas(y,X,Z=NULL,group,logodds=NULL,nfolds=10,...)
}
\arguments{
\item{y}{response variable. A vector of length p.}

\item{X}{random effect design matrix of n by p.}

\item{Z}{covariate matrix of n by q. Optional. By default, a column of 1's are added.}

\item{group}{a vector of length p, indicating the group of each variable.}

\item{logodds}{cadidate values of group level prior inclusion probability(pi) in logodds scale, a vector of length ns, default value is seq(-log10(K),log10(K),length.out = 40).}

\item{nfolds}{Number of folds being used. Default is 10.}

\item{...}{Other arguments that can be passed to \code{bivas}}

}
\details{
  \code{perf.bivas} Evaluate the performance of bivas for a given data using cross-validation.
}
\value{
The function returns a list of :
  \item{cvm}{
     The cross-validation mean of test errors.
}
  \item{cvsd}{
     The standard deviation of cross-validation test errors.
}
  \item{testErr}{
     The test errors.
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
perf.bivas(y,X,Z=NULL,group,logodds=NULL,nfolds=10,...)

}
